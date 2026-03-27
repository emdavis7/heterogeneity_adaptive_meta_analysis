#### Generate multivariate Gaussian data with required elements. ########
library(dplyr)
library(tidyr)
library(Matrix)
library(MASS) 

##################### Gaussian Response #################################
generate_gaussian_data <- function(D, p = 4, sd,
                                   data_mu = rep(0, D),
                                   corr = 0,
                                   n, beta,
                                   seed, x_seed = seed, error_seed = NA,
                                   type ="standard"){

    set.seed(x_seed)


  Gamma_tilde <- diag(kronecker(1/sd^2, rep(1,p)))

  # Give study for each row.  
  study <- c()
  for(i in 1:D){study <- append(study, rep(i, n[i]))}


  # Generate design matrix.
  if (p ==4 & type == "standard"){
    x0 <- rep(1, sum(n))
    x1 <- rnorm(sum(n), mean = 0, sd = 1)
    x2 <- rbinom(n = sum(n), size = 1, prob = .5) 
    x3 <- x1*x2
    X <- cbind(x0, x1, x2, x3)

  } else if (p ==2 & type == "standard"){
    x0 <- rep(1, sum(n))
    x1 <- rnorm(sum(n), mean = 0, sd = 1)
    X <- cbind(x0, x1)

  } else if (p > 4 & type == "standard"){
    x0 <- rep(1, sum(n))
    x1 <- rnorm(sum(n), mean = 0, sd = 1)
    x2 <- rbinom(n = sum(n), size = 1, prob = .5) 
    x3 <- x1*x2

    ## Generate remainder from a standard mv norm.
    X_remain <- mvrnorm(n=sum(n), mu=rep(0, p-4), Sigma = diag(1, nrow = (p-4)))
    X <- cbind(x0, x1, x2, x3, X_remain)

  } else if (type == "shift" & p==4){
    x1 <- vector()
    x2 <- vector()
    for(i in 1:D){
      x1 <- append(x1, rnorm(n[i], mean = 3*i, sd = 1))
      x2 <- append(x2, rbinom(n[i], size = 1, prob = (i/D - i/(2*D))))
    } 
    x3 <- x1*x2
    X <- cbind(x0, x1, x2, x3)

  } else if(type=="scale" & p==4){
    x1 <- c(rnorm(n[1], mean = 0, sd = 1)/10,
            rnorm(sum(n)-n[1], mean = 0, sd = 1))
    x2 <- rbinom(n = sum(n), size = 1, prob = .5) 
    x3 <- x1*x2
    X <- cbind(x0, x1, x2, x3)

  } else if(type == "correlated" & p==4){
    x <- mvrnorm(sum(n), mu = c(0, 0),
                 Sigma = matrix(c(1, corr, corr, 1), nrow = 2, ncol=2))
    x1 <- x[,1]
    x2 <- x[,2] 
    x3 <- x1*x2
    X <- cbind(x0, x1, x2, x3)

  }

  if (p==1){
    ## Means only.
    x0 <- rep(1, sum(n))
    X <- cbind(x0)
  }

### Create a general process for generating a single intercept and remaining covars are from mvnorm.
  if (type == "MVNorm" & p > 1){ # centered at data_mu, and then scaled up by the mu_j.
    ### Generate first column as intercept.
    x0 <- rep(1, sum(n))

    #### Generate p-1 columns from MV normal
    # create cov mat from given correlation.
    V <- corr*matrix(rep(1, (p-1)*(p-1)), nrow = (p-1),
                     ncol = (p-1)) + (1-corr)*diag(p-1)
    
    # Generate study by study:
    x_remaining_list <- lapply(1:D, function(l){
      observations <- lapply(1:n[l], function(i){
        # each observation is independent, but with corr between covariates
        c <- ifelse(data_mu[l]==0, 1, data_mu[l])
        mvrnorm(n = 1,
              mu = rep(data_mu[l], p-1), # each extra column gets same mu.
              Sigma = c^2*V) # no correlation across datasets
        })
      obs_df <- as.data.frame(do.call(rbind, observations))
    })
    x_remaining <- do.call(rbind, x_remaining_list)

    X <- as.matrix(cbind(x0, x_remaining))
  }

  ### Create block diagonal version of X
  X_list <- lapply(1:D, function(i){X[which(study==i),]})
  bdiag_X <- as.matrix(bdiag(X_list))

  ### Create block diagonal matrix of scaled XTX.
  XTX <- t(bdiag_X) %*% bdiag_X
  scaled_XTX <- Gamma_tilde %*% XTX

  ## Create a list version of scaled_XTX for alternative calculations.
  scaled_XTX_list <- lapply(1:D, function(l){
    1/sd[l]^2*
      t(X[which(study==l),])%*%
      X[which(study==l),]
  })


  #### Pre-calculate tbeta_scaled_XTX. #####
  tbeta_scaled_XTX <- t(beta) %*% scaled_XTX

  #### Generate y = XB + error
  # Add error - set seed to error_seed, if we have different ones.
  if(is.na(error_seed)==F){set.seed(error_seed)}
  error <- vector()
  for (l in 1:D){error <- append(error, rnorm(n[l], mean = 0, sd[l]))}


  y <- bdiag_X %*% beta + error

  ## Get scaled XTY.
  scaled_XTY <- Gamma_tilde %*% t(bdiag_X) %*% y

  # make the list version:
  scaled_XTY_list <- lapply(1:D, function(l){
    1/sd[l]^2*
      t(X[which(study==l),])%*%
      y[which(study==l),]
  })



  #### Pre-calculate some inputs for functions. #######################
  ## (1_d \otimes I_p) : kron
  kron <- kronecker(rep(1, D), diag(nrow = p))

  ### Also return the indv MLE:
  indv_mle <- solve(scaled_XTX) %*% scaled_XTY
  indv_mle_list <- lapply(1:D, function(l){
    solve(scaled_XTX_list[[l]]) %*% scaled_XTY_list[[l]]
  })

  ### Estimate sd to generate estimated versions of scaled xtx and xty.
  est_sd <-sapply(1:D, function(l){
    y_j <- y[which(study==l),]
    X_j <- X[which(study==l),]
    mle_j <- indv_mle_list[[l]]
    residual_j <- y_j - X_j %*% mle_j
    sse_j <- sum(residual_j^2)
    sqrt(sse_j/(n[l]-p))
  })
  
  ## Generate Gamma_tilde_est
  Gamma_tilde_est <- diag(kronecker(1/est_sd^2, rep(1,p)))
  scaled_XTX_est <- Gamma_tilde_est %*% XTX

  ## Create a list version of scaled_XTX for alternative calculations.
  scaled_XTX_est_list <- lapply(1:D, function(l){
    1/est_sd[l]^2*
      t(X[which(study==l),])%*%
      X[which(study==l),]
  })


  #### Pre-calculate tbeta_scaled_XTX. ##### 
  tbeta_scaled_XTX_est <- t(beta) %*% scaled_XTX_est



  ## Get scaled XTY.
  scaled_XTY_est <- Gamma_tilde_est %*% t(bdiag_X) %*% y

  # make the list version:
  scaled_XTY_est_list <- lapply(1:D, function(i){
    1/est_sd[i]^2*
      t(X[which(study==i),])%*%
      y[which(study==i),]
  })







  #### Use unknown SD to get combined MLE (equal pi).
  combined_mle <- solve(t(kron) %*% scaled_XTX_est %*% kron) %*% 
    t(kron) %*% scaled_XTY_est

  ## Return the optimal pi when all pi_j = pi, pi_star - sd UNKNOWN.
  var_part <-  tr( t(kron) %*% solve(scaled_XTX_est) %*% kron -
                    D * solve(t(kron) %*% scaled_XTX_est %*% kron) )
  bias_part <- t( kron %*% solve( t(kron) %*% scaled_XTX_est %*% kron) %*%
                    t(kron) %*% scaled_XTX_est %*% beta - beta) %*%
    (kron %*% solve( t(kron) %*% scaled_XTX_est %*% kron) %*% 
       t(kron) %*% scaled_XTX_est %*% beta - beta)
  pi_star <- as.numeric( var_part / (var_part + bias_part) )

  ## Return pi_star_est, using the combined and indv MLEs as plug-ins.
  bias_part_est <- t( kron %*% combined_mle - indv_mle) %*%
    (kron %*% combined_mle - indv_mle)
  pi_star_est <- as.numeric( var_part / (var_part + bias_part_est) )

  ## Use pi_star to calculate the MA KLD est
  kld_ma <- (1-pi_star)*indv_mle + pi_star*kron %*% combined_mle

  ## Use pi_star_est to calculate the MA KLD est est.
  kld_ma_est <- (1-pi_star_est)*indv_mle + pi_star_est * kron %*% combined_mle

  ## variance of the pi_star and pi_star_est estimators
  var_kld_ma_p1 <- ((1-pi_star)*diag(1, p*D) +
    pi_star * kron %*%
    solve(t(kron) %*% scaled_XTX %*% kron) %*%
    t(kron) %*% scaled_XTX)
  var_kld_ma <- var_kld_ma_p1 %*% solve(scaled_XTX) %*% t(var_kld_ma_p1)

  var_kld_ma_est_p1 <- ((1-pi_star_est)*diag(1, p*D) +
                      pi_star_est * kron %*%
                      solve(t(kron) %*% scaled_XTX %*% kron) %*%
                      t(kron) %*% scaled_XTX)
  var_kld_ma_est <- var_kld_ma_est_p1 %*% solve(scaled_XTX) %*%
    t(var_kld_ma_est_p1)


  ##########################################################
  ## Return list:
  list(study = study,
               y = y,
               bdiag_X = bdiag_X,
               X = X,
               error = error,
               scaled_XTX = scaled_XTX,
               scaled_XTY = scaled_XTY,
               scaled_XTX_list = scaled_XTX_list,
               scaled_XTY_list = scaled_XTY_list,
               tbeta_scaled_XTX = tbeta_scaled_XTX,
               # estimated SD:
              est_sd = est_sd,
              scaled_XTX_est = scaled_XTX_est,
              scaled_XTY_est = scaled_XTY_est,
              scaled_XTX_est_list = scaled_XTX_est_list,
              scaled_XTY_est_list = scaled_XTY_est_list,
              tbeta_scaled_XTX_est = tbeta_scaled_XTX_est,

               combined_mle = combined_mle,
               indv_mle = indv_mle,
               kron = kron,
               pi_star = pi_star,
               pi_star_est = pi_star_est,
               kld_ma = kld_ma,
               var_kld_ma = var_kld_ma,
               kld_ma_est = kld_ma_est,
               var_kld_ma_est = var_kld_ma_est
       )

}


