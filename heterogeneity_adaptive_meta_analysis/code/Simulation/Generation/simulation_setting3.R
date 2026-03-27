############################
## Simulation Setting 3
## Varying data generation conditions
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
library(nloptr)
library(mvmeta)

## Source functions
source(paste0(dir,"/code/functions.R"))
source(paste0(dir, "/code/generate_data.R"))



#### Conditions: ###################################
# These settings are the same across simulations.
D <- 3
p <- 4

## Create storage for settings:
n_vec <- c(20, 50, 100, 200, 500, 2000, 5000, 10000, 20000)
dgp_conditions <- 5


set.seed(3000)
beta_m <- runif(p, min = -5, max = 5)
beta_mat <- sapply(1:D, function(l){
                    round(beta_m + runif(p, min = -.25, max = .25), 2)})
beta <- c(beta_mat)

design_settings <- list( data_mus = list(c(0,0,0),
                                         c(5, 5, 5),
                                         c(10, 10, 10),
                                         c(0, 1, 2),
                                         c(0, 5, 10)),
                         sd_settings = list(c(1, 1, 1),
                                            c(5, 5, 5),
                                            c(10, 10, 10),
                                            c(1, 1, 2),
                                            c(1, 5, 10)))

###### SELECT SETTING ######

## grab the initial seed from the environment.

print(initial_seed)

results_tmp <- list()

for (design in 1:dgp_conditions){
  for (n_sel in 1:length(n_vec)){
    
  seed <- 10000*initial_seed + design*n_sel

   # get the particulars of this setting.
   n <- n_vec[n_sel]
   sd <- design_settings$sd_settings[[design]]
   data_mu <- design_settings$data_mus[[design]]


   # Run generate data
    data <- generate_gaussian_data(D = D,
                                    p = p,
                                    data_mu = data_mu,
                                    sd = sd,
                                    n = rep(n, D),
                                    beta = beta,
                                    x_seed = seed,
                                    error_seed = seed+1,
                                    type = "MVNorm")



     out <- list()
     out$n <- n
     out$design <- design
     out$data_mu <- data_mu
     out$sd <- sd
     out$beta <- beta
     out$beta_m <- beta_m
     out$seed <- seed

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


     ######### Optimization #################################################
     # For each measure of MSE,
     # (1) Use bobyqa to find minimum pi  // or, skip this step for indv/combined mle.
     # (2) Preserve convergence message
     # (3) Preserve pi
     # (4) Preserve implied beta-hat + theta-hat
     # (5) Preserve analytical MSE
     # (6) Preserve sq error (for eMSE)
     # (7) Preserve aCI UB and LB, width and coverage
     ########################################################################

     ### Measure 2: Indv MLEs ##############################################
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


     ### Measure 3: Combined MLEs ###########################################
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


     ### Measure 4: KLD Penalized with true MSE known ######################
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


     # (3a) Preserve components of c^star, for analysis.
     bPi <- blockLambda_fcn(opt_pi/max(opt_pi), p = p, D=D)

     ### Bias
     B_n <-  data$kron %*%
        solve( t(data$kron) %*% bPi %*% data$scaled_XTX %*% data$kron) %*%
           t(data$kron) %*% bPi %*% data$scaled_XTX - diag(1, p*D)
     bias_vec <- bPi %*% B_n %*% beta
     out$bias_sq <- bias_sq<- t(bias_vec) %*% bias_vec

     ### Covariance
     out$tr_cov <- tr_cov <- tr(bPi %*% B_n %*% solve(data$scaled_XTX) )

     ### Variance
     out$tr_var <- tr_var<- tr(bPi %*% B_n %*% solve(data$scaled_XTX) %*% t(B_n) %*% bPi)

     ## Implied c-star
     out$implied_cstar <- (-tr_cov)/(tr_var+bias_sq)

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


     #### Measure 5: HAM estimator #############################
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



     #### Measure 6: Biased MLE estimator #############################
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
     out$mse_biased <- mse(par = as.numeric(pi_biased),
                              scaled_XTX = data$scaled_XTX,
                              beta = beta,
                              D = D, p = p,
                              kron = data$kron)


     # (6) Sq error
     out$sqerr_biased <- as.numeric(
       t(out$beta_hat_biased - beta) %*% (out$beta_hat_biased - beta))


     # (7) biased aCI:
     biased_aCI <- kld_linear_aCI(pi = pi_biased,
                               scaled_XTX = data$scaled_XTX_est,
                               indv_mle = data$indv_mle,
                               p = p, D = D, kron = data$kron,
                               alpha = .05)
     out$biased_aCI_lb <- biased_aCI$lb
     out$biased_aCI_ub <- biased_aCI$ub

     # Record and output width and coverage:
     out$biased_aCI_width <- biased_aCI$ub - biased_aCI$lb
     out$biased_aCI_covers <- biased_aCI$lb < beta & biased_aCI$ub > beta


     ##### Measure 7: Ridgelike #####################

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



     ######## Rescaling predictors: scaled X_1,..., X_p #############
     # Calculate, for each study, the sd(X_j):
     sd_x <- lapply(1:D, function(l){
       out_vec <- c(1)
       for(i in 2:p){
         out_vec <- c(out_vec, sd(data$X[data$study==l, i]))
       }
       out_vec
     })

     # create a scaling matrix.
     inv_sd_mat <- as.matrix(bdiag(lapply(1:D,function(l){
       1/(sd_x[[l]] %*% t(sd_x[[l]]))
     })))

     # scale XTX further.
     scaled_XTX2 <- data$scaled_XTX * inv_sd_mat
     scaled_XTX_est2 <- data$scaled_XTX_est * inv_sd_mat

     # scale XTY further.
     scaled_XTY2 <- data$scaled_XTY/c(unlist(sd_x))
     scaled_XTY_est2 <- data$scaled_XTY_est/c(unlist(sd_x))

     # rescale mles.
     rescale_mle <- solve(scaled_XTX_est2) %*% scaled_XTY_est2

     # rescale beta
     rescale_beta <- beta*c(unlist(sd_x))

     ## Get the rescaled sqerr_mle.
     out$sqerr_mle_rs <- as.numeric(
       t(rescale_mle - rescale_beta) %*% (rescale_mle-rescale_beta) )

     out$mse_mle_rs <- tr(solve(scaled_XTX2))

     ### Measure 8: KLD Penalized with true MSE known, scaled #############

     # (1) use bobyqa to find the pi that minimizes MSE
     # Here, minimizing as though variance is known.
     opt2 <- bobyqa(rep(data$pi_star_est, D),
                    fn = mse,
                    scaled_XTX = scaled_XTX2,
                    beta = rescale_beta,
                    p = p, D = D, kron = data$kron,
                    lower = rep(0, D), upper = rep(1, D))

     # (2) Preserve convergence message
     out$convergence_kld_rs <- opt2$convergence

     # (3) Preserve opt_pi
     opt_pi2 <- opt2$par

     out$pi_kld1_rs <- opt_pi2[1]
     out$pi_kld2_rs <- opt_pi2[2]
     out$pi_kld3_rs <- opt_pi2[3]


     # (3a) Preserve components of c^star, for analysis.
     bPi2 <- blockLambda_fcn(opt_pi2/max(opt_pi2), p = p, D=D)

     ### Bias
     B_n2 <-  data$kron %*%
       solve( t(data$kron) %*% bPi2 %*% scaled_XTX2 %*% data$kron) %*%
       t(data$kron) %*% bPi2 %*% scaled_XTX2 - diag(1, p*D)
     bias_vec2 <- bPi2 %*% B_n2 %*% rescale_beta
     out$bias_sq_rs <- bias_sq2 <- t(bias_vec2) %*% bias_vec2

     ### Covariance
     out$tr_cov_rs <- tr_cov2 <- tr(bPi2 %*% B_n2 %*% solve(scaled_XTX2) )

     ### Variance
     out$tr_var_rs <- tr_var2 <- tr(bPi2 %*% B_n2 %*% solve(scaled_XTX2) %*%
                                      t(B_n2) %*% bPi2)

     ## Implied c-star
     out$implied_cstar_rs <- (-tr_cov2)/(tr_var2+bias_sq2)

     # (4) Preserve beta-hat and theta-hat still assuming variance known
     ests <- kld_linear_alt(pi = opt_pi2,
                            scaled_XTX = scaled_XTX2,
                            indv_mle = rescale_mle,
                            p = p, D = D, kron = data$kron)

     # Send back to original scale by dividing
     out$beta_hat_kld_rs <- ests$beta_hat #/c(unlist(sd_x))

     # (5) Preserve analytic mse
     out$mse_kld_true_rs <- mse(par = opt_pi2,
                             scaled_XTX = scaled_XTX2,
                             beta = rescale_beta, p = p,
                             D = D, kron = data$kron)

     # (6) Sq Error for eMSE
     out$sqerr_kld_rs <- as.numeric(
       t(out$beta_hat_kld_rs - rescale_beta) %*% (out$beta_hat_kld_rs - rescale_beta))



     # (7) get aCI for best pi (Sigma unknown, and using tbbeta)
     kld_aCI_rs <- kld_linear_aCI(pi = opt_pi2,
                               scaled_XTX = scaled_XTX_est2,
                               indv_mle = rescale_mle,
                               p = p, D = D, kron = data$kron,
                               alpha = .05)
     # send back to original scale:
     out$kld_aCI_lb_rs <- rs_lb <- kld_aCI_rs$lb
     out$kld_aCI_ub_rs <- rs_ub <- kld_aCI_rs$ub

     # Record and output width and coverage:
     out$kld_aCI_width_rs <- rs_ub - rs_lb
     out$kld_aCI_covers_rs <- rs_lb < rescale_beta & rs_ub > rescale_beta


     #### Measure 9: HAM estimator, with rescaling #############################
     # (1) Minimize pi for corrected UMSE (function "mse_hat")
     opt_ham2 <-  bobyqa(rep(data$pi_star_est, D),
                        fn = mse_hat,
                        scaled_XTX = scaled_XTX_est2,
                        indv_mle = rescale_mle,
                        p = p, D = D, kron = data$kron,
                        lower = rep(0, D), upper = rep(1, D))

     # (2) Preserve convergence
     out$convergence_ham_rs <- opt_ham2$convergence

     # (3) Preserve pi
     pi_ham_rs <- opt_ham2$par
     out$pi_ham1_rs <- pi_ham_rs[1]
     out$pi_ham2_rs <- pi_ham_rs[2]
     out$pi_ham3_rs <- pi_ham_rs[3]

     # (4) Preserve beta-hat and theta-hat
     ests <- kld_linear_alt(pi = pi_ham_rs,
                            scaled_XTX = scaled_XTX_est2,
                            indv_mle = rescale_mle,
                            p = p, D = D, kron = data$kron)
     out$beta_hat_ham_rs <- ests$beta_hat #/c(unlist(sd_x))

     # (5) Preserve analytical MSE.
     out$mse_ham_rs <- mse(par = as.numeric(pi_ham_rs),
                        scaled_XTX = scaled_XTX2,
                        beta = rescale_beta,
                        D = D, p = p,
                        kron = data$kron)

     # (6) Preserve sq error
     out$sqerr_ham_rs <- as.numeric(
       t(out$beta_hat_ham_rs - rescale_beta) %*% (out$beta_hat_ham_rs - rescale_beta)
     )


     # (7) HAM aCI:
     ham_aCI_rs <- kld_linear_aCI(pi = pi_ham_rs,
                               scaled_XTX = scaled_XTX_est2,
                               indv_mle = rescale_mle,
                               p = p, D = D, kron = data$kron,
                               alpha = .05)
     out$ham_aCI_lb_rs <- ham_rs_lb <- ham_aCI_rs$lb
     out$ham_aCI_ub_rs <- ham_rs_ub <- ham_aCI_rs$ub

     # Record and output width and coverage:
     out$ham_aCI_width_rs <- ham_rs_ub - ham_rs_lb
     out$ham_aCI_covers_rs <- ham_rs_lb < rescale_beta & ham_rs_ub > rescale_beta




     #### Add the parameter version of theta for the opt pi and for ham pi #####
     theta_optimal <- kld_theta(pi = opt_pi, scaled_XTX = data$scaled_XTX, beta = beta,
                                p = p, D = D, kron = data$kron)
     out$theta_opt <- theta_optimal

     theta_ham <- kld_theta(pi = pi_ham, scaled_XTX = data$scaled_XTX, beta = beta,
                            p = p, D = D, kron = data$kron)
     out$theta_ham <- theta_ham

     results_tmp <- rbind(results_tmp, out)


   }
}

results_tmp <- as.data.frame(results_tmp)

# Split into long the CI information and everything else:
ci_tmp <- results_tmp %>%
  mutate(n = as.numeric(results_tmp$n),
         design = as.character(results_tmp$design),
         sd = as.character(results_tmp$sd),
         data_mu = as.character(results_tmp$data_mu),
         seed = as.numeric(results_tmp$seed)) %>%
  dplyr::select(n, design, sd, data_mu, seed,
                mle_aCI_lb, mle_aCI_ub, mle_aCI_width, mle_aCI_covers,
                kld_aCI_lb, kld_aCI_ub, kld_aCI_width, kld_aCI_covers,
                biased_aCI_lb, biased_aCI_ub, biased_aCI_width, biased_aCI_covers,
                ham_aCI_lb, ham_aCI_ub, ham_aCI_width, ham_aCI_covers,
                rl_aCI_lb, rl_aCI_ub, rl_aCI_width, rl_aCI_covers,
                kld_aCI_lb_rs, kld_aCI_ub_rs, kld_aCI_width_rs, kld_aCI_covers_rs,
                ham_aCI_lb_rs, ham_aCI_ub_rs, ham_aCI_width_rs, ham_aCI_covers_rs,
                beta_m,  re_ma_est, theta_hat_kld, theta_hat_ham, theta_hat_biased,
                theta_opt, theta_ham, combined_mle,
                indv_mle, beta_hat_kld, beta_hat_ham, beta_hat_biased, beta_hat_rl,
                beta_hat_kld_rs, beta_hat_ham_rs,
                beta)

# For each row: transform vector to long data.frame, add parm names, and dup n and seed.
ci <- list()
for (i in 1:nrow(ci_tmp)){
  ci <- bind_rows(ci,
                  data.frame(
                    n = rep(ci_tmp$n[i], p*D),
                    sd = rep(ci_tmp$sd[i], p*D),
                    design = rep(ci_tmp$design[i],p*D),
                    data_mu = rep(ci_tmp$data_mu[i], p*D),
                    study = kronecker(1:D, rep(1, 4)),
                    parameter = rep(1:4, D),
                    seed = rep(ci_tmp$seed[i], p*D),
                    beta = unlist(ci_tmp$beta[i]),
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

                    biased_aCI_lb = unlist(ci_tmp$biased_aCI_lb[i]),
                    biased_aCI_ub = unlist(ci_tmp$biased_aCI_ub[i]),
                    biased_width = unlist(ci_tmp$biased_aCI_width[i]),
                    biased_covers = unlist(ci_tmp$biased_aCI_covers[i]),
                    rl_aCI_lb = unlist(ci_tmp$rl_aCI_lb[i]),
                    rl_aCI_ub = unlist(ci_tmp$rl_aCI_ub[i]),
                    rl_width = unlist(ci_tmp$rl_aCI_width[i]),
                    rl_covers = unlist(ci_tmp$rl_aCI_covers[i]),

                    # add rescaled;
                    kld_aCI_lb_rs = unlist(ci_tmp$kld_aCI_lb_rs[i]),
                    kld_aCI_ub_rs = unlist(ci_tmp$kld_aCI_ub_rs[i]),
                    kld_width_rs = unlist(ci_tmp$kld_aCI_width_rs[i]),
                    kld_covers_rs = unlist(ci_tmp$kld_aCI_covers_rs[i]),
                    ham_aCI_lb_rs = unlist(ci_tmp$ham_aCI_lb_rs[i]),
                    ham_aCI_ub_rs = unlist(ci_tmp$ham_aCI_ub_rs[i]),
                    ham_width_rs = unlist(ci_tmp$ham_aCI_width_rs[i]),
                    ham_covers_rs = unlist(ci_tmp$ham_aCI_covers_rs[i]),


                    # estimates and centroids
                    beta_m  = unlist(ci_tmp$beta_m[i]),
                    re_ma_est = unlist(ci_tmp$re_ma_est[i]),
                    theta_hat_kld = unlist(ci_tmp$theta_hat_kld[i]),
                    theta_hat_ham = unlist(ci_tmp$theta_hat_ham[i]),
                    theta_hat_biased = unlist(ci_tmp$theta_hat_biased[i]),
                    theta_kld = unlist(ci_tmp$theta_opt[i]),
                    theta_ham = unlist(ci_tmp$theta_ham[i]),

                    combined_mle = unlist(ci_tmp$combined_mle[i]),
                    indv_mle = unlist(ci_tmp$indv_mle[i]),
                    beta_hat_kld = unlist(ci_tmp$beta_hat_kld[i]),
                    beta_hat_ham = unlist(ci_tmp$beta_hat_ham[i]),
                    beta_hat_biased = unlist(ci_tmp$beta_hat_biased[i]),
                    beta_hat_rl = unlist(ci_tmp$beta_hat_rl[i]),
                    beta_hat_kld_rs = unlist(ci_tmp$beta_hat_kld_rs[i]),
                    beta_hat_ham_rs = unlist(ci_tmp$beta_hat_ham_rs[i])

                  ))
}
# We will want to output this ci dataframe.


results_no_ci <- results_tmp %>% dplyr::select(-c(mle_aCI_lb, mle_aCI_ub, mle_aCI_width, mle_aCI_covers,
                                                    kld_aCI_lb, kld_aCI_ub, kld_aCI_width, kld_aCI_covers,
                                                    biased_aCI_lb, biased_aCI_ub, biased_aCI_width, biased_aCI_covers,
                                                    ham_aCI_lb, ham_aCI_ub, ham_aCI_width, ham_aCI_covers,
                                                    rl_aCI_lb, rl_aCI_ub, rl_aCI_width, rl_aCI_covers,
                                                  kld_aCI_lb_rs, kld_aCI_ub_rs, kld_aCI_width_rs, kld_aCI_covers_rs,
                                                  ham_aCI_lb_rs, ham_aCI_ub_rs, ham_aCI_width_rs, ham_aCI_covers_rs,
                                                    beta_m, re_ma_est, theta_hat_kld, theta_hat_ham, theta_hat_biased, theta_opt, theta_ham, combined_mle,
                                                    indv_mle, beta_hat_kld, beta_hat_ham, beta_hat_biased, beta_hat_rl,
                                                  beta_hat_kld_rs, beta_hat_ham_rs,
                                                    #kld_scaled_bias, ham_scaled_bias,
                                                    beta))

results_tmp2 <- as.data.frame(apply(
  dplyr::select(results_no_ci, -c(data_mu, sd)), 2, as.numeric))

results_tmp2$data_mu <- as.character(results_no_ci$data_mu)
results_tmp2$sd <- as.character(results_no_ci$sd)

# Calculate as a ratio of MLE's MSE.
results_tmp2$mse_ratio_rl <- results_tmp2$mse_rl/results_tmp2$mse_mle

results_tmp2$mse_ratio_kld <- results_tmp2$mse_kld_true/results_tmp2$mse_mle
results_tmp2$mse_ratio_ham <- results_tmp2$mse_ham /results_tmp2$mse_mle
results_tmp2$mse_ratio_kld_rs <- results_tmp2$mse_kld_true_rs/results_tmp2$mse_mle_rs
results_tmp2$mse_ratio_ham_rs <- results_tmp2$mse_ham_rs /results_tmp2$mse_mle_rs

results_tmp2$mse_ratio_biased <- results_tmp2$mse_biased /results_tmp2$mse_mle

# Compare to the scenario where we assume total homogeneity
results_tmp2$mse_ratio2_rl <- results_tmp2$mse_rl/results_tmp2$mse_combined

results_tmp2$mse_ratio2_kld <- results_tmp2$mse_kld_true/results_tmp2$mse_combined
results_tmp2$mse_ratio2_ham<- results_tmp2$mse_ham/results_tmp2$mse_combined
results_tmp2$mse_ratio2_biased <- results_tmp2$mse_biased /results_tmp2$mse_combined

### Write the generated results to RDS
if(home_comp_flag==F){
saveRDS(results_tmp2,
        file = paste0(dir,
                      "/simulation_data/setting3_mse/results_",
                      initial_seed,".rds"),
        compress = TRUE)
saveRDS(ci,
        file = paste0(dir,
                      "/simulation_data/setting3_ci/results_",
                      initial_seed,".rds"),
        compress = TRUE)
}
