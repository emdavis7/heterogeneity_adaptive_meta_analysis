############################
# Simulation: Setting 4
### Condition 2: 5 outlier studies with Beta_1 from Unif(-.25, .25)
############################

library(dplyr)
library(tidyr)
library(Matrix)
library(mvmeta)
library(nloptr)
library(extraDistr)

## Source functions
source(paste0(dir,"/code/functions.R"))
source(paste0(dir, "/code/generate_data.R"))


### Fixed conditions ####
D <- 20
sd <- 1
p <- 5
q_j <- rep(5, D)
n_j <- rep(200, D)


# Generate beta_{1j}
set.seed(9876)
beta_1_common <- runif(1, min = -5, max = 5)
beta_1_1 <- sapply(1:5, function(l){
  round(runif(1, min = -.25, max = .25), 2)
})
beta_1_2 <- sapply(6:D, function(l){
  round(beta_1_common + runif(1, min = -.25, max = .25), 2)
})
beta_1 <- c(beta_1_1, beta_1_2)

# Generate remaining betas, and assemble beta
beta_list <- lapply(1:D, function(l){
  q <- q_j[l]
  beta_remaining <- rnorm(q-1)
  beta_1_j <- beta_1[l]
  append(beta_remaining, beta_1_j, after = 1)
})
beta <- unlist(beta_list)


out <- list()
results_tmp <- list()


################## Begin Replicates: ##########################
# Roll through replicates:
for (i in 1:1000){
  if(i%%100==0){print(i)}
  set.seed(i)
  out$seed <- i
  out$beta_1 <- beta_1

  data <- generate_gaussian_data(D = D,
                                 p = p,
                                 sd = rep(sd, D),
                                 n = n_j,
                                 beta = beta,
                                 x_seed = i,
                                 error_seed = i + 1,
                                 type = "MVNorm")


  ### Isolate data of interest. ####
  # indv_mle[2] for each dataset
  beta1_mle <- matrix(data$indv_mle, nrow = p)[2,]

  subcov_XTX_est_list <- lapply(1:D, function(l){
    data$scaled_XTX_est_list[[l]][2,2]
  })
  subcov_XTZ_est_list <- lapply(1:D, function(l){
    data$scaled_XTX_est_list[[l]][2,-2]
  })
  subcov_ZTX_est_list <- lapply(1:D, function(l){
    data$scaled_XTX_est_list[[l]][-2,2]
  })
  subcov_ZTZ_est_list <- lapply(1:D, function(l){
    data$scaled_XTX_est_list[[l]][-2,-2]
  })


  ### Create scaled X_j'M_j X_j #####
  scaled_XTMX_est_list <- lapply(1:D, function(l){
    subcov_XTX_est_list[[l]] -
      subcov_XTZ_est_list[[l]] %*% solve(subcov_ZTZ_est_list[[l]]) %*% subcov_ZTX_est_list[[l]]
  })


  # Block diagonal:
  scaled_XTMX_est <- diag(scaled_XTMX_est_list)
  beta1_var <- solve(scaled_XTMX_est)

  # Not estimated
  scaled_XTMX_list <- lapply(1:D, function(l){
    data$est_sd[l]*(subcov_XTX_est_list[[l]] -
                      subcov_XTZ_est_list[[l]] %*% solve(subcov_ZTZ_est_list[[l]]) %*% subcov_ZTX_est_list[[l]])
  })
  scaled_XTMX <- diag(scaled_XTMX_list)

  ######### Optimization #################################################
  # For each measure of MSE,
  # (1) Use bobyqa to find minimum pi  // or, skip this step for indv/combined mle.
  # (2) Preserve convergence message
  # (3) Preserve pi
  # (4) Preserve implied beta-hat (no theta-hat for this sim study)
  # (5) Preserve analytical MSE
  # (6) Preserve sq error (for eMSE)
  # (7) Preserve aCI UB and LB, width and coverage
  ########################################################################

  ### Measure 0: RE MA, just on single var ####
  re_ma <- mvmeta(beta1_mle, S = diag(beta1_var))

  # (4) Preserve estimate
  out$re_ma_est <- re_ma$coefficients[1]

  # (6) Preserve eMSE
  out$sqerr_rema <- as.numeric(
    t(out$re_ma_est - beta_1) %*% (out$re_ma_est - beta_1))

  # Sq error against common:
  out$sqerr_rema_common <- as.numeric(
    t(out$re_ma_est - beta_1_common) %*% (out$re_ma_est - beta_1_common))

  # (7) Preserve confidence intervals.
  re_ma_lb <- re_ma$coefficients[1] - 1.96*sqrt(re_ma$vcov[1,1])
  re_ma_ub <- re_ma$coefficients[1] + 1.96*sqrt(re_ma$vcov[1,1])

  re_ma_width <- re_ma_ub - re_ma_lb
  re_ma_covers <- beta_1 >= re_ma_lb & beta_1 <= re_ma_ub


  out$re_ma_lb <- re_ma_lb
  out$re_ma_ub <- re_ma_ub
  out$re_ma_width <- re_ma_ub - re_ma_lb
  out$re_ma_covers_common <- beta_1_common >= re_ma_lb & beta_1_common <= re_ma_ub
  out$re_ma_covers <- beta_1 >= re_ma_lb & beta_1 <= re_ma_ub

  ### Measure 1: Indv MLEs #####
  # (4) preserve MLE
  out$indv_mle <- beta1_mle

  # (5) analytical mse
  out$mse_mle <- tr(beta1_var)

  # (6) sqerr mle
  out$sqerr_mle <- as.numeric(
    t(beta1_mle - beta_1) %*% (beta1_mle - beta_1) )

  out$mle_aCI_lb <- sapply(1:length(beta1_mle), function(l){
    beta1_mle[l] - qnorm(.975, mean = 0, sd = 1)*sqrt(diag(beta1_var)[l])
  })

  out$mle_aCI_ub <- sapply(1:length(beta1_mle), function(l){
    beta1_mle[l] + qnorm(.975, mean = 0, sd = 1)*sqrt(diag(beta1_var)[l])
  })

  # Record and output width and coverage:
  out$mle_aCI_width <- out$mle_aCI_ub - out$mle_aCI_lb
  out$mle_aCI_covers <- out$mle_aCI_lb < beta_1 & out$mle_aCI_ub > beta_1

  ### Measure 2: KLD Optimization, true MSE ####
  # Find pi to minimize MSE
  # Use .2 as starting value: in situation with unmatched covariates, we may not be able to calculate pi_star_est.
  kron <- kronecker(rep(1, D), diag(nrow = 1))
  opt1 <- bobyqa(rep(.2, D),
                 fn = mse,
                 scaled_XTX = scaled_XTMX,
                 beta = beta_1,
                 p = 1, D = D, kron = kron,
                 lower = rep(0, D),
                 upper = rep(1, D))

  # (2) Preserve convergence message
  out$convergence_kld <- opt1$convergence
  # (3) Preserve opt_pi
  out$pi_kld <-  opt_pi <- opt1$par

  # (3a) Check trVar, trCov, sq_bias:
  bPi <- blockLambda_fcn(opt_pi/max(opt_pi), p = 1, D=D)
  stack <- rep(1, D)
  ### Bias
  B_n <-  stack %*%
    solve( t(stack) %*% bPi %*% scaled_XTMX %*% stack ) %*%
    t(stack) %*% bPi %*% scaled_XTMX - diag(1, 1*D)
  bias_vec <- bPi %*% B_n %*% beta_1
  out$bias_sq <- bias_sq<- t(bias_vec) %*% bias_vec

  ### Covariance
  out$tr_cov <- tr_cov <- tr(bPi %*% B_n %*% solve(scaled_XTMX) )

  ### Variance
  out$tr_var <- tr_var<- tr(bPi %*% B_n %*% solve(scaled_XTMX) %*% t(B_n) %*% bPi)

  ## Implied c-star
  out$implied_cstar <- (-tr_cov)/(tr_var+bias_sq)

  # (4) Preserve beta-hat and theta-hat still assuming variance known
  est_kld <-  kld_linear_alt(pi = opt1$par,
                             scaled_XTX = scaled_XTMX,
                             indv_mle = beta1_mle,
                             p = 1, D = D, kron = kron)
  out$beta_hat_kld <- est_kld$beta_hat
  out$theta_hat_kld <- est_kld$theta_hat

  # (5) Preserve analytic mse
  out$mse_kld_true <- opt1$value

  # (6) Sq Error for eMSE
  # Sq Error for eMSE
  out$sqerr_kld <- as.numeric(
    t(est_kld$beta_hat - beta_1) %*% (est_kld$beta_hat - beta_1))

  # Sq error against common:
  out$sqerr_kld_common <- as.numeric(
    t(est_kld$theta_hat - beta_1_common) %*% (est_kld$theta_hat - beta_1_common))


  # (7) get aCI for best pi (Sigma unknown, and using tbbeta)
  kld_aCI <- kld_linear_aCI(pi = opt1$par,
                            scaled_XTX = scaled_XTMX_est,
                            indv_mle = beta1_mle,
                            p = 1, D = D, kron = kron,
                            alpha = .05)

  out$kld_aCI_lb <- kld_aCI$lb
  out$kld_aCI_ub <- kld_aCI$ub

  # Record and output width and coverage:
  out$kld_aCI_width <- kld_aCI$ub - kld_aCI$lb
  out$kld_aCI_covers <- kld_aCI$lb < beta_1 & kld_aCI$ub > beta_1

  ##### Measure 3: HAM ###############
  # Find pi to minimize pseudo MSE.
  pi_ham_opt <-  bobyqa(rep(.2, D),
                        fn = mse_hat,
                        scaled_XTX = scaled_XTMX_est,
                        indv_mle = beta1_mle,
                        p = 1, D = D, kron = kron,
                        lower = rep(0, D),
                        upper = rep(1, D) )

  # (2) Preserve convergence
  out$convergence_ham <- pi_ham_opt$convergence

  # (3) Preserve pi
  out$pi_ham <- pi_ham_opt$par

  # (4) Preserve beta-hat and theta-hat
  est_ham <- kld_linear_alt(pi = pi_ham_opt$par,
                            scaled_XTX = scaled_XTMX_est,
                            indv_mle = beta1_mle,
                            p = 1, D = D, kron = kron)
  out$beta_hat_ham <- est_ham$beta_hat
  out$theta_hat_ham <- est_ham$theta_hat

  # (5) Preserve analytical MSE.
  out$mse_ham <- mse(par = pi_ham_opt$par,
                     scaled_XTX = scaled_XTMX,
                     beta = beta_1,
                     p = 1, D = D, kron = kron)

  # (6) Preserve sq error
  # Sq Error for eMSE
  out$sqerr_ham <- as.numeric(
    t(est_ham$beta_hat - beta_1) %*% (est_ham$beta_hat - beta_1))

  # Sq error against common:
  out$sqerr_ham_common <- as.numeric(
    t(est_ham$theta_hat - beta_1_common) %*% (est_ham$theta_hat - beta_1_common))


  # (7) HAM aCI:
  ham_aCI <- kld_linear_aCI(pi = pi_ham_opt$par,
                            scaled_XTX = scaled_XTMX_est,
                            indv_mle = beta1_mle,
                            p = 1, D = D, kron = kron,
                            alpha = .05)

  out$ham_aCI_lb <- ham_aCI$lb
  out$ham_aCI_ub <- ham_aCI$ub

  # Record and output width and coverage:
  out$ham_aCI_width <- ham_aCI$ub - ham_aCI$lb
  out$ham_aCI_covers <- ham_aCI$lb < beta_1 & ham_aCI$ub > beta_1



  #### Measure 4: Ridge-like ####
  ## (1) Choose lambda for ridgelike to minimize umse.
  opt_rl <- bobyqa(1,
                   fn = rl_umse,
                   scaled_XTX = scaled_XTMX_est,
                   indv_mle = beta1_mle,
                   p = 1, D = D,
                   lower = 0)

  # (2) Preserve convergence
  out$convergence_rl <- opt_rl$convergence

  # (3) Preserve lambda
  out$lambda_rl <- opt_rl$par

  # (4) Preserve beta-hat
  out$beta_hat_rl <- ridgelike(lambda = out$lambda_rl,
                               scaled_XTX = scaled_XTMX_est,
                               scaled_XTY = scaled_XTMX_est %*% beta1_mle,
                               p = 1, D = D)

  # (5) Analytical MSE for RL
  out$mse_rl <- rl_mse(lambda = out$lambda_rl,
                       scaled_XTX = scaled_XTMX,
                       beta = beta_1,
                       p = 1, D = D)

  # (6) sq error for rl
  out$sqerr_rl <- as.numeric(
    t(out$beta_hat_rl - beta_1) %*% (out$beta_hat_rl - beta_1) )


  # (7) RL aCI:
  rl_aCI_out <- rl_aCI(lambda = out$lambda_rl,
                       scaled_XTX = scaled_XTMX_est,
                       indv_mle = beta1_mle,
                       p = 1, D = D,
                       alpha = .05)
  out$rl_aCI_lb <- rl_aCI_out$lb
  out$rl_aCI_ub <- rl_aCI_out$ub

  # Record and output width and coverage:
  out$rl_aCI_width <- rl_aCI_out$ub - rl_aCI_out$lb
  out$rl_aCI_covers <- rl_aCI_out$lb < beta_1 & rl_aCI_out$ub > beta_1

  results_tmp <- rbind(results_tmp, out)
}

#### Tabulate results #####
results_tmp <- as.data.frame(results_tmp)

# Split into long the CI information and everything else:
ci_tmp <- results_tmp %>%
  mutate(seed = as.numeric(results_tmp$seed)) %>%
  dplyr::select(seed,
                re_ma_covers,
                mle_aCI_lb, mle_aCI_ub, mle_aCI_width, mle_aCI_covers,
                kld_aCI_lb, kld_aCI_ub, kld_aCI_width, kld_aCI_covers,
                ham_aCI_lb, ham_aCI_ub, ham_aCI_width, ham_aCI_covers,
                rl_aCI_lb, rl_aCI_ub, rl_aCI_width, rl_aCI_covers,
                indv_mle, beta_hat_kld, beta_hat_ham, beta_hat_rl,
                beta_1,
                pi_kld, pi_ham)
# For each row: transform vector to long data.frame, add parm names, and dup seed.
ci <- list()
for (i in 1:nrow(ci_tmp)){
  ci <- bind_rows(ci,
                  data.frame(
                    seed = rep(ci_tmp$seed[i], D),
                    re_ma_covers = unlist(ci_tmp$re_ma_covers[i]),
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
                    indv_mle = unlist(ci_tmp$indv_mle[i]),
                    beta_hat_kld = unlist(ci_tmp$beta_hat_kld[i]),
                    beta_hat_ham = unlist(ci_tmp$beta_hat_ham[i]),
                    beta_hat_rl = unlist(ci_tmp$beta_hat_rl[i]),
                    beta_1 = unlist(ci_tmp$beta_1[i]),

                    # pi
                    pi_kld = unlist(ci_tmp$pi_kld[i]),
                    pi_ham = unlist(ci_tmp$pi_ham[i]),

                    study = 1:D)
  )
}
# We will want to output this ci dataframe.

results_no_ci <- results_tmp %>% dplyr::select(-c(re_ma_covers,
                                                  kld_aCI_lb, kld_aCI_ub, kld_aCI_width, kld_aCI_covers,
                                                  ham_aCI_lb, ham_aCI_ub, ham_aCI_width, ham_aCI_covers,
                                                  rl_aCI_lb, rl_aCI_ub, rl_aCI_width, rl_aCI_covers,
                                                  mle_aCI_lb, mle_aCI_ub, mle_aCI_width, mle_aCI_covers,
                                                  indv_mle, beta_hat_kld, beta_hat_ham, beta_hat_rl,
                                                  pi_kld, pi_ham,
                                                  beta_1))

results_tmp2 <- as.data.frame(apply(results_no_ci, 2, as.numeric))

# Calculate as a ratio of MLE's MSE.
results_tmp2$mse_ratio_opt_kld <- results_tmp2$mse_kld_true/results_tmp2$mse_mle
results_tmp2$mse_ratio_ham <- results_tmp2$mse_ham/results_tmp2$mse_mle
results_tmp2$mse_ratio_rl <- results_tmp2$mse_rl/results_tmp2$mse_mle


## Save to simulation data folder.
saveRDS(results_tmp2, file = "./data/results_setting4_rep_level2.RDS")
saveRDS(ci, file = "./data/results_setting4_coeff_level2.RDS")

