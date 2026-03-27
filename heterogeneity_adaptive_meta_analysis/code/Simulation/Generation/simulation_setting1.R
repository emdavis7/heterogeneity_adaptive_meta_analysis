############################
# Simulation: Setting 1
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
library(MASS)
library(mvmeta)
library(nloptr)
library(Hmisc)

## Source functions
source(paste0(dir,"/code/functions.R"))
source(paste0(dir, "/code/generate_data.R"))




#### Conditions: ###################################
set.seed(1)

# These settings are the same across simulations.
D <- 3
p_vec <- c(2, 4, 10, 20)
sd <- rep(1, D)


# Create list of n combinations.
n <- t(matrix(c(100, 100, 100,
                100, 200, 100,
                100, 300, 100,
                100, 200, 200,
                100, 300, 200,
                100, 300, 300,
                200, 200, 200,
                200, 300, 200,
                200, 300, 300,
                300, 300, 300
                ), ncol = 10))

# gen settings list to select while looping.
sel_df <- expand_grid(p_settings = 1:length(p_vec),
                      n_settings = 1:nrow(n))

results_tmp <- list()

for (set in 1:nrow(sel_df)){
   ## using setting to dictate the random seed.
   seed <- 10000*initial_seed + set

   # get the particulars of this setting.
   select_n <- sel_df$n_settings[set]
   n_set <- as.numeric(n[select_n,])

   select_p <- sel_df$p_settings[set]
   p <- as.numeric(p_vec[select_p])

   set.seed(987654)
   beta_m <- runif(p, min = -5, max = 5)
   beta_mat <- sapply(1:D, function(l){round(beta_m + runif(p, min = -.25, max = .25), 2)})
   beta <- c(beta_mat)

   out <- list()


   # Run generate data
   data <- generate_gaussian_data(D = D,
                                  p = p,
                                  sd = sd,
                                  n = n_set,
                                  beta = beta,
                                  x_seed = seed,
                                  error_seed = seed+1)

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


   out$re_ma_est <- as.vector(re_ma$coefficients)
   out$sqerr_re_ma <-  as.numeric(
      t(data$kron %*% out$re_ma_est - beta) %*% (data$kron %*% out$re_ma_est - beta)
   )
   ####################################################


   #### Alternative heterogeneity measures #####
   studywise_eucdist <- sapply(1:D, function(l){
     diff <- beta_mat[,l] - beta_m
     sqrt(t(diff) %*% diff)
   })
   out$eucdist1 <- studywise_eucdist[1]
   out$eucdist2 <- studywise_eucdist[2]
   out$eucdist3 <- studywise_eucdist[3]

   #######################################################

   out$n_set <- paste0(n_set, collapse = ", ")
   out$p <- p
   out$beta <- beta
   out$beta_m <- beta_m
   out$set <- set
   out$seed <- seed


   ######### Optimization #################################################
   # For each measure of MSE,
   # (1) Use bobyqa to find minimum pi  // skip this step for indv/combined mle.
   # (2) Preserve convergence message
   # (3) Preserve pi
   # (4) Preserve implied beta-hat + theta-hat
   # (5) Preserve analytical MSE
   # (6) Preserve sq error (for eMSE)
   # (7) Preserve aCI UB and LB, width and coverage
   # (8) Preserve t adj aCI UB and LB, width and coverage.
   ########################################################################

   ### Measure 1: Indv MLEs ##############################################
   # (4) implied beta
   out$indv_mle <- data$indv_mle

   # (5) Analytical MSE for MLE
   out$mse_mle <- tr(solve(data$scaled_XTX))

   # (6) Sum of squared error for eMSE
   out$sqerr_mle <- as.numeric(
      t(data$indv_mle - beta) %*% (data$indv_mle-beta) )

   # (7) For comparison, aCI
   out$mle_aCI_lb <- data$indv_mle -
     qnorm(.975)*sqrt(diag(solve(data$scaled_XTX_est)))
   out$mle_aCI_ub <- data$indv_mle +
     qnorm(.975)*sqrt(diag(solve(data$scaled_XTX_est)))

   # Record and output width and coverage:
   out$mle_aCI_width <- out$mle_aCI_ub - out$mle_aCI_lb
   out$mle_aCI_covers <- out$mle_aCI_lb < beta & out$mle_aCI_ub > beta

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



   # (7) get aCI for best pi (Sigma unknown, and using tbbeta)
   kld_aCI <- kld_linear_aCI(pi = opt_pi,
                             scaled_XTX = data$scaled_XTX_est,
                             indv_mle = data$indv_mle,
                             p = p, D = D, kron = data$kron,
                             alpha = .05)
   out$kld_aCI_lb <- kld_aCI$lb
   out$kld_aCI_ub <- kld_aCI$ub

   # Record and output width and coverage:
   out$kld_aCI_width <- kld_aCI$ub - kld_aCI$lb
   out$kld_aCI_covers <- kld_aCI$lb < beta & kld_aCI$ub > beta


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


   # (7) HAM aCI:
   ham_aCI <- kld_linear_aCI(pi = pi_ham,
                              scaled_XTX = data$scaled_XTX_est,
                              indv_mle = data$indv_mle,
                              p = p, D = D, kron = data$kron,
                              alpha = .05)
   out$ham_aCI_lb <- ham_aCI$lb
   out$ham_aCI_ub <- ham_aCI$ub

   # Record and output width and coverage:
   out$ham_aCI_width <- ham_aCI$ub - ham_aCI$lb
   out$ham_aCI_covers <- ham_aCI$lb < beta & ham_aCI$ub > beta


   ##### Measure 5: Ridgelike #####################

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
   out$beta_hat_rl <- ridgelike(lambda = out$lambda_rl,
                          scaled_XTX = data$scaled_XTX_est,
                          scaled_XTY = data$scaled_XTY_est,
                          p = p, D = D)

   # (5) Analytical MSE for RL
   out$mse_rl <- rl_mse(lambda = out$lambda_rl,
                                scaled_XTX = data$scaled_XTX,
                                beta = beta,
                                p = p, D = D)

   # (6) sq error for rl
   out$sqerr_rl <- as.numeric(
      t(out$beta_hat_rl - beta) %*% (out$beta_hat_rl - beta) )


   # (7) RL aCI:
   rl_aCI_out <- rl_aCI(lambda = out$lambda_rl,
                        scaled_XTX = data$scaled_XTX_est,
                        indv_mle = data$indv_mle,
                        p = p, D = D,
                        alpha = .05)
   out$rl_aCI_lb <- rl_aCI_out$lb
   out$rl_aCI_ub <- rl_aCI_out$ub

   # Record and output width and coverage:
   out$rl_aCI_width <- rl_aCI_out$ub - rl_aCI_out$lb
   out$rl_aCI_covers <- rl_aCI_out$lb < beta & rl_aCI_out$ub > beta




   #### Add the parameter version of theta for the opt pi and for ham pi #####
   theta_optimal <- kld_theta(pi = opt_pi, scaled_XTX = data$scaled_XTX, beta = beta,
                              p = p, D = D, kron = data$kron)
   out$theta_opt <- theta_optimal

   theta_ham <- kld_theta(pi = pi_ham, scaled_XTX = data$scaled_XTX, beta = beta,
                              p = p, D = D, kron = data$kron)
   out$theta_ham <- theta_ham

   results_tmp <- rbind(results_tmp, out)


}


results_tmp <- as.data.frame(results_tmp)

# Split into long the CI information and everything else:
ci_tmp <- results_tmp %>%
  mutate(n = as.character(results_tmp$n_set),
         p = as.numeric(results_tmp$p),
         seed = as.numeric(results_tmp$seed)) %>%
  dplyr::select(n, p, seed,
                mle_aCI_lb, mle_aCI_ub, mle_aCI_width, mle_aCI_covers,
                kld_aCI_lb, kld_aCI_ub, kld_aCI_width, kld_aCI_covers,
                ham_aCI_lb, ham_aCI_ub, ham_aCI_width, ham_aCI_covers,
                rl_aCI_lb, rl_aCI_ub, rl_aCI_width, rl_aCI_covers,
                beta_m, re_ma_est, theta_hat_kld, theta_hat_ham, theta_opt, theta_ham, combined_mle,
                indv_mle, beta_hat_kld, beta_hat_ham, beta_hat_rl,
                beta)

# For each row: transform vector to long data.frame, add parm names, and dup n and seed.
ci <- list()
for (i in 1:nrow(ci_tmp)){
  ci <- bind_rows(ci,
                  data.frame(
                   n = rep(ci_tmp$n[i], length(unlist(ci_tmp$kld_aCI_lb[i]))),
                   p = rep(ci_tmp$p[i], length(unlist(ci_tmp$kld_aCI_lb[i]))),
                   seed = rep(ci_tmp$seed[i], length(unlist(ci_tmp$kld_aCI_lb[i]))),
                   mle_aCI_lb = unlist(ci_tmp$mle_aCI_lb[i]),
                   mle_aCI_ub = unlist(ci_tmp$mle_aCI_ub[i]),
                   mle_width = unlist(ci_tmp$mle_aCI_width[i]),
                   mle_covers = unlist(ci_tmp$mle_aCI_covers[i]),
                   kld_aCI_lb = unlist(ci_tmp$kld_aCI_lb[i]),
                   kld_aCI_ub = unlist(ci_tmp$kld_aCI_ub[i]),
                   kld_width = unlist(ci_tmp$kld_aCI_width[i]),
                   kld_covers = unlist(ci_tmp$kld_aCI_covers[i]),
                   ham_aCI_lb = unlist(ci_tmp$ham_aCI_lb[i]),
                   ham_aCI_ub = unlist(ci_tmp$ham_aCI_ub[i]),
                   ham_width = unlist(ci_tmp$ham_aCI_width[i]),
                   ham_covers = unlist(ci_tmp$ham_aCI_covers[i]),
                   rl_aCI_lb = unlist(ci_tmp$rl_aCI_lb[i]),
                   rl_aCI_ub = unlist(ci_tmp$rl_aCI_ub[i]),
                   rl_width = unlist(ci_tmp$rl_aCI_width[i]),
                   rl_covers = unlist(ci_tmp$rl_aCI_covers[i]),


                   # estimates and centroids
                   beta_m  = unlist(ci_tmp$beta_m[i]),
                   re_ma_est = unlist(ci_tmp$re_ma_est[i]),
                   theta_hat_kld = unlist(ci_tmp$theta_hat_kld[i]),
                   theta_hat_ham = unlist(ci_tmp$theta_hat_ham[i]),
                   theta_kld = unlist(ci_tmp$theta_opt[i]),
                   theta_ham = unlist(ci_tmp$theta_ham[i]),
                   combined_mle = unlist(ci_tmp$combined_mle[i]),
                   indv_mle = unlist(ci_tmp$indv_mle[i]),
                   beta_hat_kld = unlist(ci_tmp$beta_hat_kld[i]),
                   beta_hat_ham = unlist(ci_tmp$beta_hat_ham[i]),
                   beta_hat_rl = unlist(ci_tmp$beta_hat_rl[i]),
                   beta = unlist(ci_tmp$beta[i]),

                   study = kronecker(1:D, rep(1, ci_tmp$p[i])),
                   parameter = rep(1:ci_tmp$p[i], D)
             ))
}
# We will want to output this ci dataframe.

results_no_ci <- results_tmp %>% dplyr::select(-c(kld_aCI_lb, kld_aCI_ub, kld_aCI_width, kld_aCI_covers,
                                                  ham_aCI_lb, ham_aCI_ub, ham_aCI_width, ham_aCI_covers,
                                                  rl_aCI_lb, rl_aCI_ub, rl_aCI_width, rl_aCI_covers,
                                                  mle_aCI_lb, mle_aCI_ub, mle_aCI_width, mle_aCI_covers,
                                                  beta_m, re_ma_est, theta_hat_kld, theta_hat_ham, theta_opt, theta_ham, combined_mle,
                                                  indv_mle, beta_hat_kld, beta_hat_ham, beta_hat_rl,
                                                  beta))

results_tmp2 <- as.data.frame(apply(
  dplyr::select(results_no_ci, -n_set), 2, as.numeric))

results_tmp2$n <- as.character(results_no_ci$n_set)

# Calculate as a ratio of MLE's MSE.
results_tmp2$mse_ratio_rl <- results_tmp2$mse_rl/results_tmp2$mse_mle
results_tmp2$mse_ratio_kld <- results_tmp2$mse_kld_true/results_tmp2$mse_mle
results_tmp2$mse_ratio_ham <- results_tmp2$mse_ham /results_tmp2$mse_mle

# Compare to the scenario where we assume total homogeneity
results_tmp2$mse_ratio2_rl <- results_tmp2$mse_rl/results_tmp2$mse_combined
results_tmp2$mse_ratio2_kld <- results_tmp2$mse_kld_true/results_tmp2$mse_combined
results_tmp2$mse_ratio2_ham<- results_tmp2$mse_ham/results_tmp2$mse_combined


### Write the generated results to RDS
if(home_comp_flag==F){
saveRDS(results_tmp2,
        file = paste0(dir,
                      "/simulation_data/setting1_mse/results_",
                      initial_seed,".rds"),
        compress = TRUE)
saveRDS(ci,
        file = paste0(dir,
                      "/simulation_data/setting1_ci/results_",
                      initial_seed,".rds"),
        compress = TRUE)
}
