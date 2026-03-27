############################
# Simulation: Setting 2 
############################

# Switch depending on computer:
home_comp_flag <- T

if (home_comp_flag==T){
  dir <- getwd()
  # nrep * number of jobs = number of sim reps.
  nrep <- 1
  initial_seed <- 1
} else {
  dir <-  "/SET_DIRECTORY_HERE"
  initial_seed <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
  nrep <- 1
}


library(dplyr)
library(tidyr)
library(Matrix)
library(mvmeta)
library(nloptr)

## Source functions
source(paste0(dir,"/code/functions.R"))
source(paste0(dir, "/code/generate_data.R"))




#### Conditions: ###################################
# These settings are the same across simulations.
p <- 4
n <- 200
sd <- 1

## Create storage for settings:
conditions <- 4
D_list <- c(rep(5, conditions),
            rep(10, conditions),
            rep(15, conditions))
cond_list <- rep(c(1:conditions), 3)

betas <- c("beta_mat <- sapply(1:D, function(l){signif(beta_m, 2)})
                     beta <- c(beta_mat)",
           "beta_mat <- sapply(1:D, function(l){
                       signif(mvrnorm(1, beta_m, .1*diag(1, p)), 2)})
                     beta <- c(beta_mat)",
           "beta_mat <- sapply(1:D, function(l){
                       signif(mvrnorm(1, beta_m, .5*diag(1, p)), 2)})
                     beta <- c(beta_mat)",
           "beta_mat <- sapply(1:D, function(l){
                       if(l < 4){signif(beta_m, 2)
                         } else {signif(mvrnorm(1, beta_m, diag(1, p)), 2)}
                     })
                     beta <- c(beta_mat)"
           )

results_tmp <- list()
out <- list()

#### Initialize vectors for long output: ########
studywise_eucdist <- vector()
d_to_pair <- vector()
cond_to_pair  <- vector()
seed_to_pair <- vector()
study_ham_mse <- vector()
study_mle_mse <- vector()

beta_to_pair <- vector()
beta_hat_kld <- vector()
kld_aCI_lb <- vector()
kld_aCI_ub  <- vector()
kld_aCI_width <- vector()
kld_aCI_covers <- vector()
beta_hat_ham <- vector()
ham_aCI_lb <- vector()
ham_aCI_ub  <- vector()
ham_aCI_width <- vector()
ham_aCI_covers <- vector()
indv_mle <- vector()
re_ma_est <- vector()
mle_aCI_lb <- vector()
mle_aCI_ub  <- vector()
mle_aCI_width <- vector()
mle_aCI_covers <- vector()
combined_mle <- vector()
beta_hat_rl <- vector()


##### Rolls through settings. ########
for (set in 1:length(cond_list)){

 ## using setting to dictate the random seed.
 seed <- 10000*initial_seed + set

 # get the particulars of this setting.
 D <- as.numeric(D_list[set])

 # Set the shared center:
 set.seed(101)
 beta_m <- signif(runif(p, min = -5, max = 5), 2)

 # Randomly select new draws of individual study parms.
 set.seed(seed)
 beta_cond <- cond_list[set]
 eval(parse(text = betas[beta_cond]))
 bTb <- t(beta) %*% beta

 # Run generate data
 data <- generate_gaussian_data(D = D,
                                p = p,
                                sd = rep(sd, D),
                                n = rep(n, D),
                                beta = beta,
                                x_seed = seed,
                                error_seed = seed+1)



  out$D <- D
  out$condition  <- cond_list[set]
  out$seed <- seed
  out$set <- set

   ################ Estimate MV I^2 ###################
   mle_mat <- t(matrix(data$indv_mle, nrow = p))
   cov_list <- lapply(1:D, function(l){solve(data$scaled_XTX_est_list[[l]])})
   C_F <- solve(Reduce(`+`, data$scaled_XTX_est_list))

   re_ma <- mvmeta(mle_mat, S = cov_list)
   C_R <- re_ma$vcov
   I_sq <- ( det(C_R)^(1/p) - det(C_F)^(1/p) )/det(C_R)^(1/p)
   Q_pvalue <- qtest(re_ma)$pvalue[1]

   out$I_sq <- I_sq
   out$Q_pvalue <- Q_pvalue

   re_ma_tmp <- as.vector(re_ma$coefficients)
   re_ma_est <- c(re_ma_est, re_ma_tmp)
   out$sqerr_re_ma <-  as.numeric(
      t(data$kron %*% re_ma_tmp - beta) %*% (data$kron %*% re_ma_tmp - beta)
   )
   ####################################################


   ####Studywise Euclidean Distance#####
   studywise_eucdist <- c(studywise_eucdist, sapply(1:D, function(l){
     diff <- beta_mat[,l] - beta_m
     sqrt(t(diff) %*% diff)
   }) )

   d_to_pair <- c(d_to_pair, 1:D)
   cond_to_pair <- c(cond_to_pair, rep(cond_list[set], D) )
   seed_to_pair <- c(seed_to_pair, rep(seed, D))
   beta_to_pair <- c(beta_to_pair, beta)

   ######### Optimization #################################################
   # For each measure of MSE,
   # (1) Use bobyqa to find minimum pi  // or, skip this step for indv/combined mle.
   # (2) Preserve convergence message
   # (3) Preserve pi
   # (4) Preserve implied beta-hat (no theta-hat for this sim study)
   # (5) Preserve analytical MSE
   # (5b) Study-specific aMSE
   # (6) Preserve sq error (for eMSE)
   # (7) Preserve aCI UB and LB, width and coverage
   ########################################################################

   ### Measure 1: Indv MLEs #####
   #  (4) preserve beta
   indv_mle <- c(indv_mle, data$indv_mle)

   # (5) Preserve analytical MSE
   out$mse_mle <- tr(solve(data$scaled_XTX))

   # (5b) Get study MLE MSE:
   study_mle_mse <- c(study_mle_mse,
                      sapply(1:D, function(l){
                        tr(solve(data$scaled_XTX_list[[l]]))
                      })
   )

   # (6) sum of sq err for eMSE
   out$sqerr_mle <- as.numeric(
     t(data$indv_mle - beta) %*% (data$indv_mle-beta) )

   # (7) Preserve aCI UB and LB, width and coverage
   mle_lb_tmp <- sapply(1:length(data$indv_mle), function(l){
     data$indv_mle[l] - qnorm(.975, mean = 0, sd = 1)*sqrt(diag(solve(data$scaled_XTX_est))[l])
   })
   mle_aCI_lb <- c(mle_aCI_lb, mle_lb_tmp)
   mle_ub_tmp <- sapply(1:length(data$indv_mle), function(l){
     data$indv_mle[l] + qnorm(.975, mean = 0, sd = 1)*sqrt(diag(solve(data$scaled_XTX_est))[l])
   })
   mle_aCI_ub  <- c(mle_aCI_ub, mle_ub_tmp)
   mle_width_tmp <- mle_ub_tmp - mle_lb_tmp
   mle_aCI_width <- c(mle_aCI_width, mle_width_tmp)
   mle_covers_tmp <- mle_ub_tmp > beta & mle_lb_tmp < beta
   mle_aCI_covers <- c(mle_aCI_covers, mle_covers_tmp)


   ### Measure 2: Combined MLE ####
   # (4) Implied beta vector
   combined_mle <- c(combined_mle, data$kron %*% data$combined_mle)

   # (5) Analytical MSE for combined MLE
   combined_bias <- (data$kron %*% solve(t(data$kron) %*% data$scaled_XTX %*% data$kron) %*%
                       t(data$kron) %*% data$scaled_XTX - diag(1, nrow = p*D))%*%beta
   combined_var <-  solve(t(data$kron) %*% data$scaled_XTX %*% data$kron)
   out$mse_combined <- t(combined_bias) %*% combined_bias + tr(combined_var)

   # (6) Sum of squared error, for eMSE
   out$sqerr_combined_MLE <- as.numeric(
     t(data$kron %*% data$combined_mle - beta) %*% (data$kron %*% data$combined_mle - beta)
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
   out$opt_pi <- opt_pi <- opt1$par

   # (4) Preserve beta-hat still assuming variance known
   ests <- kld_linear_alt(pi = opt_pi,
                          scaled_XTX = data$scaled_XTX,
                          indv_mle = data$indv_mle,
                          p = p, D = D, kron = data$kron)
   beta_hat_kld <- c(beta_hat_kld, ests$beta_hat)
   # kld_theta_hat <- c(kld_theta_hat, ests$theta_hat)

   # (5) Preserve analytic mse
   out$mse_kld_true <- mse(par = opt_pi,
                           scaled_XTX = data$scaled_XTX,
                           beta = beta, p = p,
                           D = D, kron = data$kron)

   # (6) Sq Error for eMSE
   out$sqerr_kld <- as.numeric(
     t(ests$beta_hat - beta) %*% (ests$beta_hat - beta))



   # (7) get aCI for best pi (Sigma unknown, and using tbbeta)
   kld_aCI <- kld_linear_aCI(pi = opt1$par,
                             scaled_XTX = data$scaled_XTX_est,
                             indv_mle = data$indv_mle,
                             p = p, D = D, kron = data$kron,
                             alpha = .05)
   kld_aCI_lb <- c(kld_aCI_lb, kld_aCI$lb)
   kld_aCI_ub <- c(kld_aCI_ub, kld_aCI$ub)
   kld_aCI_width <- c(kld_aCI_width, kld_aCI$ub - kld_aCI$lb)
   kld_aCI_covers <- c(kld_aCI_covers, kld_aCI$lb < beta & kld_aCI$ub > beta)


   #### Measure 4: HAM estimator #############################
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
   out$pi_ham <- pi_ham <- opt_ham$par


   # (4) Preserve beta-hat and theta-hat
   ests <- kld_linear_alt(pi = pi_ham,
                          scaled_XTX = data$scaled_XTX_est,
                          indv_mle = data$indv_mle,
                          p = p, D = D, kron = data$kron)
   beta_hat_ham <- c(beta_hat_ham, ests$beta_hat)
   # out$theta_hat_ham <- ests$theta_hat

   # (5) Preserve analytical MSE.
   out$mse_ham <- mse(par = as.numeric(pi_ham),
                      scaled_XTX = data$scaled_XTX,
                      beta = beta,
                      D = D, p = p,
                      kron = data$kron)

    # (5b) Get study MSEs:
    study_ham_mse <- c(study_ham_mse,
                       studywise_mse(pi = as.numeric(pi_ham),
                                     scaled_XTX = data$scaled_XTX,
                                     scaled_XTX_list = data$scaled_XTX_list,
                                     beta = beta,
                                     beta_mat = beta_mat,
                                     p = p, D = D, kron = data$kron) )

    # (6) Preserve sq error
    out$sqerr_ham <- as.numeric(
      t(ests$beta_hat - beta) %*% (ests$beta_hat - beta))


    # (7) HAM aCI
    ham_aCI <- kld_linear_aCI(pi = pi_ham,
                             scaled_XTX = data$scaled_XTX_est,
                             indv_mle = data$indv_mle,
                             p = p, D = D, kron = data$kron,
                             alpha = .05)
    ham_aCI_lb <- c(ham_aCI_lb, ham_aCI$lb)
    ham_aCI_ub <- c(ham_aCI_ub, ham_aCI$ub)
    ham_aCI_width <- c(ham_aCI_width, ham_aCI$ub - ham_aCI$lb)
    ham_aCI_covers <- c(ham_aCI_covers, ham_aCI$lb < beta & ham_aCI$ub > beta)


   ######## Measure 5: Ridgelike ##########################

    ## (1) Choose lambda for ridgelike to minimize umse.
    opt_rl <- bobyqa(1,
                     fn = rl_umse,
                     scaled_XTX = data$scaled_XTX_est,
                     indv_mle = data$indv_mle,
                     p = p, D = D,
                     lower = 0)

    # (2) Preserve convergence
    out$convergence_rl <- opt_rl$convergence

    # (3) Preserve lambda
    out$lambda_rl <- opt_rl$par

    # (4) Preserve beta-hat
    rl_est <- ridgelike(lambda = out$lambda_rl,
                        scaled_XTX = data$scaled_XTX_est,
                        scaled_XTY = data$scaled_XTY_est,
                        p = p, D = D)
    beta_hat_rl <- c(beta_hat_rl, rl_est)

    # (5) Analytical MSE for RL
    out$mse_rl <- rl_mse(lambda = out$lambda_rl,
                         scaled_XTX = data$scaled_XTX,
                         beta = beta,
                         p = p, D = D)

    # (6) sq error for rl
    out$sqerr_rl <- as.numeric(
      t(rl_est - beta) %*% (rl_est - beta) )



   results_tmp <- rbind(results_tmp, out)


 }

##### Replicate level output: #########
results_tmp <- as.data.frame(results_tmp)

##### study-level output dataset: #######
study_level <- data.frame(studywise_eucdist, D = d_to_pair,
                      beta_cond = cond_to_pair, seed = seed_to_pair,
                      study_mle_mse,
                      study_ham_mse)


##### Coeff-level output: ######
coeff_level <- data.frame(D = kronecker(d_to_pair, rep(1, 4) ),
                          beta_cond = kronecker(cond_to_pair, rep(1,4)),
                          seed = kronecker(seed_to_pair, rep(1,4)),
                          parameter = 1:4,
                          beta = beta_to_pair,
                          beta_hat_kld,
                          kld_aCI_lb ,
                          kld_aCI_ub  ,
                          kld_aCI_width ,
                          kld_aCI_covers ,
                          beta_hat_ham,
                          ham_aCI_lb ,
                          ham_aCI_ub  ,
                          ham_aCI_width ,
                          ham_aCI_covers,
                          indv_mle,
                          re_ma_est,
                          mle_aCI_lb ,
                          mle_aCI_ub  ,
                          mle_aCI_width ,
                          mle_aCI_covers,
                          combined_mle)

results_tmp2 <- as.data.frame(apply(dplyr::select(results_tmp, -c(opt_pi, pi_ham)), 2, as.numeric))
results_tmp2$opt_pi <- results_tmp$opt_pi
results_tmp2$pi_ham <- results_tmp$pi_ham

# Calculate as a ratio of MLE's MSE.
results_tmp2$mse_ratio_rl <- results_tmp2$mse_rl/results_tmp2$mse_mle
results_tmp2$mse_ratio_kld <- results_tmp2$mse_kld_true/results_tmp2$mse_mle
results_tmp2$mse_ratio_ham <- results_tmp2$mse_ham/results_tmp2$mse_mle

# Calculate as a ratio of MLE's MSE.
results_tmp2$mse_ratio2_rl <- results_tmp2$mse_rl/results_tmp2$mse_combined
results_tmp2$mse_ratio2_kld <- results_tmp2$mse_kld_true/results_tmp2$mse_combined
results_tmp2$mse_ratio2_ham <- results_tmp2$mse_ham/results_tmp2$mse_combined


### Write the generated results to RDS
if(home_comp_flag==F){
saveRDS(results_tmp2,
        file = paste0(dir,
                      "/simulation_data/setting2_rep_level/results_",
                      initial_seed,".rds"),
        compress = TRUE)
saveRDS(study_level,
        file = paste0(dir,
                      "/simulation_data/setting2_study_level/results_",
                      initial_seed,".rds"),
        compress = TRUE)
saveRDS(coeff_level,
        file = paste0(dir,
                      "/simulation_data/setting2_coeff_level/results_",
                      initial_seed,".rds"),
        compress = TRUE)
}
