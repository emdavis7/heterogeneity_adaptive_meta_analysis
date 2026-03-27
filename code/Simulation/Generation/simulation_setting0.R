#############header###############
# Simulation: Setting 0
#############header end###############

# Switch depending on computer:
home_comp_flag <- T

if (home_comp_flag==T){
  dir <- getwd() 
  nrep <- 200
 initial_seed <- 1
} else {
  dir <-  "/SET_DIRECTORY_HERE"
  initial_seed <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
  nrep <- 1
}

library(dplyr)
library(tidyr)
library(Matrix)
library(MASS)
library(nloptr)
library(Hmisc)

## Source functions
source(paste0(dir,"/code/functions.R"))
source(paste0(dir, "/code/generate_data.R"))



#### Simulation Conditions: ###################################
set.seed(1)

# These settings are the same across simulations.
D <- 3
p <- 4
sd <- rep(.5, D)

set.seed(987654)
beta_m <- runif(p, min = -5, max = 5)
beta_mat <- sapply(1:D, function(l){signif(mvrnorm(1, beta_m, .1*diag(1, p)), 2)})
beta <- c(beta_mat)

# Create list of n combinations.
n <- t(matrix(c(100, 100, 100,
                100, 500, 100,
                100, 500, 500,
                500, 500, 500), ncol = 4))


#### Begin loop over settings ########
results_tmp <- list()
out <- list()

for (set in 1:nrow(n)){
  loop_seed <- 100000*initial_seed + set

  # get the particulars of this setting.
  n_set <- as.numeric(n[set,])

  # Seed for generating the same X.
  external_seed <- 1
   # Putting this as set ensures same X is generated across
   # machines when sample sizes are equal.

    # Run generate data
    data <- generate_gaussian_data(D = D,
                                    p = p,
                                    sd = sd,
                                    n = n_set,
                                    beta = beta,
                                    x_seed = external_seed,
                                    error_seed = loop_seed)


     out$n_set <- paste0(n_set, collapse = ", ")
     out$external_seed <- external_seed
     out$seed <- loop_seed
     out$beta <- beta
     out$beta_m <- beta_m

 ######### Contents for Each Measure##############################################
 # For each measure of MSE,
 # (1) Use bobyqa to find minimum pi  // or, skip this step for indv/combined mle.
 # (2) Preserve convergence message
 # (3) Preserve pi
 # (4) Preserve implied beta-hat + theta-hat
 # (5) Preserve analytical MSE
 # (6) Preserve sq error (for eMSE)
 #################end note#######################################################

     ### Measure 1: Indv MLEs ##############################################
     # (4) implied beta
     out$indv_mle <- data$indv_mle

     # (5) Analytical MSE for MLE
     out$mse_mle <- tr(solve(data$scaled_XTX))

     # (6) Sum of squared error for eMSE
     out$sqerr_mle <- as.numeric(
      t(data$indv_mle - beta) %*% (data$indv_mle-beta) )

     ### Measure 2: Combined MLEs ###########################################
     # (4) Implied beta vector
     out$combined_mle <- data$kron %*% data$combined_mle

     # (5) Analytical MSE for combined MLE
     combined_bias <- (data$kron %*% solve(t(data$kron) %*% data$scaled_XTX %*% data$kron) %*%
                         t(data$kron) %*% data$scaled_XTX - diag(1, nrow = p*D))%*%beta
     combined_var <-  solve(t(data$kron) %*% data$scaled_XTX %*% data$kron)
     out$mse_combined <- t(combined_bias) %*% combined_bias + tr(combined_var)

     # (6) Sum of squared error, for eMSE
     out$sqerr_combined_MLE <- as.numeric(
       t(out$combined_mle - beta) %*% (out$combined_mle - beta)
     )

     ### Measure 3: KLD Penalized with true MSE known ######################
     # (1) use bobyqa to find the pi that minimizes MSE
     # Here, minimizing as though variance is known.
     opt1 <- bobyqa(rep(data$pi_star_est, D),
                    fn = mse,
                    scaled_XTX = data$scaled_XTX,
                    beta = beta,
                    p = p, D = D, kron = data$kron,
                    lower = rep(0, D), upper = rep(1, D))

     # (2) Preserve convergence message
     out$convergence_kld <- opt1$convergence

     # (3) Preserve opt_pi
     opt_pi <- opt1$par

     out$pi_kld1 <- opt_pi[1]
     out$pi_kld2 <- opt_pi[2]
     out$pi_kld3 <- opt_pi[3]

     # (4) Preserve beta-hat and theta-hat still assuming variance known
     ests <- kld_linear_alt(pi = opt_pi,
                    scaled_XTX = data$scaled_XTX,
                    indv_mle = data$indv_mle,
                    p = p, D = D, kron = data$kron)
     out$beta_hat_kld <- ests$beta_hat
     out$theta_hat_kld <- ests$theta_hat

     # (5) Preserve analytic mse
     out$mse_kld_true <- mse(par = opt_pi,
                        scaled_XTX = data$scaled_XTX,
                        beta = beta, p = p,
                        D = D, kron = data$kron)

     # (6) Sq Error for eMSE
     out$sqerr_kld <- as.numeric(
       t(out$beta_hat_kld - beta) %*% (out$beta_hat_kld - beta))

     ### Measure 4: KLD penalization, minimizing UMSE #########################
     # (1) Find pi to minimize UMSE
     umse_opt <- bobyqa(rep(data$pi_star_est, D),
                        fn = sure_mse,
                        scaled_XTX = data$scaled_XTX_est,
                        indv_mle= data$indv_mle,
                        p = p, D = D, kron = data$kron,
                        lower = rep(0, D), upper = rep(1, D))

     # (2) preserve convergence message
     out$convergence_umse <- umse_opt$convergence

     # (3) Preserve pi UMSE
     pi_umse <- umse_opt$par

     out$pi_umse1 <- pi_umse[1]
     out$pi_umse2 <- pi_umse[2]
     out$pi_umse3 <- pi_umse[3]


     # (4) Get beta-hat (variance now is unknown)
     ests <- kld_linear_alt(pi = pi_umse,
                           scaled_XTX = data$scaled_XTX_est,
                           indv_mle = data$indv_mle,
                           p = p, D = D, kron = data$kron)
     out$beta_hat_umse <- ests$beta_hat
     out$theta_hat_umse <- ests$theta_hat

     # (5) Analytic MSE for the pi that minimizes UMSE
     out$mse_at_min_umse <- mse(par = umse_opt$par,
                                 scaled_XTX = data$scaled_XTX,
                                 beta = beta,
                                 p = p, D = D,
                                 kron = data$kron)

     # (6) Sq error
     out$sqerr_umse <- as.numeric(
       t(out$beta_hat_umse - beta) %*% (out$beta_hat_umse - beta))



     #### Measure 5: Minimize the biased MSE ################################
     # (1) Find pi
     opt_biased <- bobyqa(rep(data$pi_star_est, D),
                            fn = naive_mse,
                            scaled_XTX = data$scaled_XTX_est,
                            indv_mle = data$indv_mle,
                            p = p, D = D, kron = data$kron,
                            lower = rep(0, D), upper = rep(1, D))

     # (2) Preserve convergence
     out$convergence_biased <- opt_biased$convergence

     # (3) Preserve pi
     pi_biased <- opt_biased$par
     out$pi_biased1 <- pi_biased[1]
     out$pi_biased2 <- pi_biased[2]
     out$pi_biased3 <- pi_biased[3]

     # (4) Preserve beta-hat and theta-hat at this pi vector;
     ests <- kld_linear_alt(pi = pi_biased,
                            scaled_XTX = data$scaled_XTX_est,
                            indv_mle = data$indv_mle,
                            p = p, D = D, kron = data$kron)
     out$beta_hat_biased <- ests$beta_hat
     out$theta_hat_biased <- ests$theta_hat

     # (5) Analytical mse at this minimum:
     out$mse_at_biased <- mse(par = as.numeric(pi_biased),
                                scaled_XTX = data$scaled_XTX,
                                beta = beta,
                                D = D, p = p,
                                kron = data$kron)

     # (6) Sq error
     out$sqerr_biased <- as.numeric(
       t(out$beta_hat_biased - beta) %*% (out$beta_hat_biased - beta))

     #### Measure 6: HAM estimator #############################
     # (1) Minimize pi for corrected UMSE (function "mse_hat")
     opt_ham <-  bobyqa(rep(data$pi_star_est, D),
                                 fn = mse_hat,
                                 scaled_XTX = data$scaled_XTX_est,
                                 indv_mle = data$indv_mle,
                                 p = p, D = D, kron = data$kron,
                                 lower = rep(0, D), upper = rep(1, D))

     # (2) Preserve convergence
     out$convergence_ham <- opt_ham$convergence

     # (3) Preserve pi
     pi_ham <- opt_ham$par
     out$pi_ham1 <- pi_ham[1]
     out$pi_ham2 <- pi_ham[2]
     out$pi_ham3 <- pi_ham[3]

     # (4) Preserve beta-hat and theta-hat
     ests <- kld_linear_alt(pi = pi_ham,
                            scaled_XTX = data$scaled_XTX_est,
                            indv_mle = data$indv_mle,
                            p = p, D = D, kron = data$kron)
     out$beta_hat_ham <- ests$beta_hat
     out$theta_hat_ham <- ests$theta_hat

     # (5) Preserve analytical MSE.
     out$mse_ham <- mse(par = as.numeric(pi_ham),
                                 scaled_XTX = data$scaled_XTX,
                                 beta = beta,
                                 D = D, p = p,
                                 kron = data$kron)

     # (6) Preserve sq error
     out$sqerr_ham <- as.numeric(
       t(out$beta_hat_ham - beta) %*% (out$beta_hat_ham - beta)
     )


     results_tmp <- rbind(results_tmp, out)
}


#### Reorganize results and add ratios #######
results_tmp <- as.data.frame(results_tmp)
results_tmp2 <- results_tmp %>% mutate_if(all.is.numeric, as.numeric)
results_tmp2$n <- as.character(results_tmp$n_set)

# Calculate as a ratio of MLE's MSE.
results_tmp2$mse_ratio_opt_kld <- results_tmp2$mse_kld_true/results_tmp2$mse_mle
results_tmp2$mse_ratio_umse <- results_tmp2$mse_at_min_umse/results_tmp2$mse_mle
results_tmp2$mse_ratio_biased <- results_tmp2$mse_at_biased/results_tmp2$mse_mle
results_tmp2$mse_ratio_ham <- results_tmp2$mse_ham/results_tmp2$mse_mle


# Calculate as a ratio of the combined MLE's MSE:
results_tmp2$mse_ratio2_opt_kld <- results_tmp2$mse_kld_true/results_tmp2$mse_combined
results_tmp2$mse_ratio2_umse <- results_tmp2$mse_at_min_umse/results_tmp2$mse_combined
results_tmp2$mse_ratio2_biased <- results_tmp2$mse_at_biased/results_tmp2$mse_combined
results_tmp2$mse_ratio2_ham <- results_tmp2$mse_ham/results_tmp2$mse_combined

### Write the generated results to RDS #######
if(home_comp_flag==F){
saveRDS(results_tmp2,
        file = paste0(dir,
                      "/simulation_data/sim0/results_",
                      initial_seed,".rds"),
        compress = TRUE)
}
