#######
# Functions for Simulation
#######
library(Matrix)
library(dplyr)
library(numDeriv)
library(purrr)

################################### Auxiliary Functions ##########################################3

### Create auxiliary elements that will help with computation. #########

### Trace: ####
tr <- function(x){
  sum(diag(x))
}


### Combined MSE ####
mse_mle_combined <- function(scaled_XTX, beta, p, D, kron){
  A <-  solve(t(kron) %*% scaled_XTX %*% kron) %*%
    t(kron) %*% scaled_XTX
  combined_bias <- (kron %*% A - diag(1, nrow = p*D)) %*% beta
  combined_var <- solve(t(kron) %*% scaled_XTX %*% kron)
  t(combined_bias) %*% combined_bias + D*tr(combined_var)
}



########################################################
## Functions to create e_i and e_i e_i'
## e_i e_i' : selBlock, a D times D matrix with 1 in the (i,i) position and 0 elsewhere
## e_i : e_i, D times 1 with 1 in ith pos and 0 elsewhere.
########################################################

## e_i e_i' : selBlock, a function returning D times D matrix with 1 in the (i,i) position and 0 elsewhere
selBlock <- function(i, D){
  replace(numeric(D), i, 1) %*% t(replace(numeric(D), i, 1))
}

## e_i : e_i, D times 1 with 1 in ith pos and 0 elsewhere.
e_i <- function(i, D){
  replace(numeric(D), i, 1)
}


### Estimate sigma_sq for each dataset. #####
# If sigma_sq is unknown, estimate for each study using sum of square errors.
est_sigma_sq <- function(data, p){
  x_matrix <- data$X
  y <- data$y
  residual <- (diag(nrow(x_matrix)) - x_matrix %*%
                 solve(t(x_matrix) %*% x_matrix) %*%
                 t(x_matrix)) %*% y
  residual_sq <- residual^2
  colnames(residual_sq) <- "residual_sq"
  tmp <- data.frame(residual_sq, study = data$study)
  tmp2 <- tmp %>% group_by(study) %>% summarize(sse = sum(residual_sq),
                                                adjust = n() - p)
  return(tmp2$sse/tmp2$adjust)
}

####### Lambda Block ####
# Takes a D by 1 lambda vector and turns it into a block diagonal matrix pD by pD
blockLambda_fcn <- function(lambda, p=4, D){
    # Diagonalize - uses hadamard product to keep only the diagonal.
    Lambda <- lambda %*% matrix(data = rep(1, D), nrow = 1, ncol = D) * diag(D)
    kronecker(Lambda, diag(p))
}


###### True MSE #######
mse <- function(par, scaled_XTX, beta, D, p, kron, par_type = "pi"){
  if (par_type == "lambda") {
    pi <- par/(par+1)
  } else {
    pi <- par
  }

  # Block diagonal pi.
  bPi <- blockLambda_fcn(pi, p = p, D = D)


  if (all.equal(pi, rep(0, D))==T){tr(solve(scaled_XTX))
  } else{

    # Theta maker:
    A <- solve( t(kron) %*% bPi %*% scaled_XTX %*% kron) %*%
      t(kron) %*% bPi %*% scaled_XTX
    # bias_comp
    b_comp <- bPi %*% (kron %*% A - diag(1, p*D))
    bias_sq <- t(beta) %*% t(b_comp) %*% b_comp %*% beta

    var_mle <- tr(solve(scaled_XTX))
    var_betahat <- tr( t(b_comp) %*% b_comp %*% solve(scaled_XTX))
    cov_piece <- 2*tr(
      solve(scaled_XTX) %*% t(b_comp)
    )

    # Return the sum.
    var_mle +
      cov_piece +
      var_betahat +
      bias_sq
  }
}


###### SURE MSE #######
# Parameterizing with pi for now.
sure_mse <- function(par, scaled_XTX, indv_mle, D, p, kron, par_type = "pi"){
  if (par_type == "lambda") {
    pi <- par/(par+1)
  } else {
    pi <- par
  }

  # Catch all zero:
  if (all.equal(pi, rep(0, D))==T){tr(solve(scaled_XTX))
  } else{

  # Block diagonal pi.
  bPi <- blockLambda_fcn(pi, p = p, D = D)

  # Theta maker:
  A <- solve( t(kron) %*% bPi %*% scaled_XTX %*% kron) %*%
    t(kron) %*% bPi %*% scaled_XTX
  # bias_comp
  b_comp <- bPi %*% (kron %*% A - diag(1, p*D))

  var_mle <- tr(solve(scaled_XTX))

  diff_sq <- t(indv_mle) %*% t(b_comp) %*% b_comp %*% indv_mle

  cov_piece <- 2*tr(
    solve(scaled_XTX) %*% t(b_comp)
  )

  # Return the sum.
  var_mle +
    diff_sq +
    cov_piece
  }
}


### MSE Hat without C:
mse_hat <- function(par, scaled_XTX, indv_mle, D, p, kron, par_type = "pi", c_comp = NULL){
  if (par_type == "lambda") {
    pi <- par/(par+1)
  } else {
    pi <- par
  }

  # Catch all zero:
  if (all.equal(pi, rep(0, D))==T){tr(solve(scaled_XTX))
  } else{
    # Reparameterize to include scaling:
    c <- max(pi)
    pi2 <- pi/c

    # Block diagonal pi.
    bPi <- blockLambda_fcn(pi2, p = p, D = D)

    # Theta maker:
    A <- solve( t(kron) %*% bPi %*% scaled_XTX %*% kron) %*%
      t(kron) %*% bPi %*% scaled_XTX
    # bias_comp
    b_comp <- bPi %*% (kron %*% A - diag(1, p*D))

    var_mle <- tr(solve(scaled_XTX))

    diff_sq <- t(indv_mle) %*% t(b_comp) %*% b_comp %*% indv_mle

    var_betahat <- tr( t(b_comp) %*% b_comp %*% solve(scaled_XTX))

    cov_piece <- tr(
      solve(scaled_XTX) %*% t(b_comp)
    )

    naive_cstar <- -cov_piece /(diff_sq + var_betahat)
    hc_star <- -cov_piece/diff_sq

    if (is.null(c_comp)){c_comp <-  naive_cstar}# If c_comp is blank, use c_naive.

       c^2*diff_sq -
       2*c*diff_sq*(naive_cstar)
  }
}




# Naive MSE (BMSE):
naive_mse <- function(par, scaled_XTX, indv_mle, D, p, kron, par_type = "pi"){
  if (par_type == "lambda") {
    pi <- par/(par+1)
  } else {
    pi <- par
  }

  # Block diagonal pi.
  bPi <- blockLambda_fcn(pi, p = p, D = D)

  # Theta maker:
  A <- solve( t(kron) %*% bPi %*% scaled_XTX %*% kron) %*%
    t(kron) %*% bPi %*% scaled_XTX

  var_mle <- tr(solve(scaled_XTX))

  diff_sq <- t(indv_mle - kron %*% A %*% indv_mle) %*%
    bPi %*% bPi %*%
    (indv_mle - kron %*% A %*% indv_mle)

  cov_piece <- 2*tr(
    solve(scaled_XTX) %*% (- bPi + bPi %*% kron %*% A
    )
  )

  extra_piece1 <- bPi %*%  (kron %*%A - diag(1, p*D))
  extra_piece <- tr(extra_piece1 %*% solve(scaled_XTX) %*% t(extra_piece1))

  # Return the sum.
  var_mle + diff_sq + cov_piece + extra_piece

}






########### Dataset-wise Functions ##########

studywise_mse <- function(pi, scaled_XTX, scaled_XTX_list, beta, beta_mat, p, D, kron){
  # Block diagonal pi.
  bPi <- blockLambda_fcn(pi, p = p, D = D)
  
  if (all.equal(pi, rep(0, D))==T){
    sapply(1:D, function(l){
      tr(solve(scaled_XTX_list[[l]]))
    })
  } else{
    # Theta maker:
    A <- solve( t(kron) %*% bPi %*% scaled_XTX %*% kron) %*%
      t(kron) %*% bPi %*% scaled_XTX
    
    # Returns a vector of MSEs.
    sapply(1:D, function(j){
      S_j <- kronecker(selBlock(j, D), diag(1, p))
      pi_j <- pi[j]
      scaled_XTX_j <- scaled_XTX_list[[j]]
      beta_j <- beta_mat[, j]
      var_comp <- (1-pi_j)^2 * solve(scaled_XTX_j) +
        pi_j^2 * A %*% solve(scaled_XTX) %*% t(A) +
        2*(1-pi_j)*pi_j* t(kron) %*% S_j %*% solve(scaled_XTX) %*% t(A)
      bias_comp <- pi_j*(A %*% beta - beta_j)
      
      # Return the sum of trace of var and sq bias.
      tr(var_comp) + t(bias_comp) %*% bias_comp
    })
  }
}




###### cstar ###############
c_star <- function(const_bpi, scaled_XTX, p, D, beta, kron){

  # Turn constrained pi vec into block diag matrix;
  bpi <- blockLambda_fcn(const_bpi, p = p, D = D)
  K <- kron
  inv_scaled_XTX <- solve(scaled_XTX)

  # big I, pk times pk;
  big_I  <- diag(1, nrow= p*D, ncol = p*D)

  # Theta Maker;
  A <- solve( t(K) %*% bpi %*% scaled_XTX %*% K ) %*% t(K) %*% bpi %*% scaled_XTX
  B <- K %*% A - big_I

  # trCov;
  trCov <- tr( inv_scaled_XTX %*% t(B) %*% bpi)

  # trVar;
  trVar <- tr(bpi %*% B %*% inv_scaled_XTX %*% t(B) %*% bpi)

  # sqBias;
  sqBias <- t(beta) %*% t(B) %*% t(bpi) %*% bpi %*% B %*% beta

  min( -trCov / (trVar + sqBias), 1 )
}


########## Functions for Ridgelike #######################

## Create the contrast matrix ####
contr_matrix <- function(D, p){
  choices <- c(1, -1)

  # Then, every combination of 1 - d
  main <- combn(c(1:D), 2)

  # Create a d by ncol(main) matrix of zeros.
  mat <- matrix( rep(0, D*ncol(main)), nrow = D, ncol = ncol(main))

  # Then, the first row of main gives us the position of 1, and
  # the second gives position of negative 1.
  ones <- main[1,]
  neg <- main[2,]
  for (j in 1:ncol(main)){
    mat[ones[j], j] <- 1
    mat[neg[j], j] <- -1
  }
  kronecker(mat %*% t(mat), diag(1, p))
}


### Ridge-like Estimator #####
ridgelike <- function(lambda, scaled_XTX, scaled_XTY,
                        p, D){
  cont <- contr_matrix(D, p)
  solve(scaled_XTX + 2*lambda*cont) %*% scaled_XTY
}


rl_mse <- function(lambda, scaled_XTX, p, D, beta){
  cont <- contr_matrix(D, p)
  interior <- solve(scaled_XTX + 2*lambda*cont)
  var_comp <- tr(interior %*% scaled_XTX %*% interior)

  bias <- (interior %*% scaled_XTX - diag(1, p*D)) %*% beta
  bias_comp <- t(bias) %*% bias
  var_comp + bias_comp
}

rl_umse <- function(lambda, scaled_XTX, p, D, indv_mle){
  cont <- contr_matrix(D, p)
  interior <- solve(scaled_XTX + 2*lambda*cont)
  R <- interior %*% scaled_XTX

  bias <- (R - diag(1, nrow = p*D)) %*% indv_mle
  bias_sq <- t(bias) %*% bias

  var_comp <- tr((2*R - diag(1, nrow= p*D))*solve(scaled_XTX))
  bias_sq + var_comp
}


rl_aCI <- function(lambda, scaled_XTX, p, D, indv_mle, alpha = .05){
  cont <- contr_matrix(D, p)
  interior <- solve(scaled_XTX + 2*lambda*cont)
  R <- interior %*% scaled_XTX

  beta_hat <- R %*% indv_mle
  var <- R %*% solve(scaled_XTX) %*% t(R)

  z_crit <- qnorm(1-alpha/2, mean = 0, sd = 1)
  lb <- beta_hat - z_crit * sqrt(diag(var))
  ub <- beta_hat + z_crit * sqrt(diag(var))
  list(lb = lb, ub = ub)
}


################ KLD linear #############################################################


## Alternative expression for linear KLD:
kld_linear_alt <- function(pi, scaled_XTY = NULL, scaled_XTX, indv_mle, p, D, kron){

  # Make block diagonal pi;
  bPi <- blockLambda_fcn(pi, p, D)

  
  
  theta_hat <- solve(t(kron) %*% bPi %*%
          scaled_XTX %*% kron) %*%
    t(kron) %*% bPi %*%
    scaled_XTX %*% indv_mle

  beta_hat <- (diag(1, nrow = p*D)-bPi) %*% indv_mle +
    bPi %*% kron %*%
    theta_hat


  list(beta_hat = beta_hat, theta_hat = theta_hat)
}


kld_theta <- function(pi, scaled_XTX, beta, p, D, kron){

  # Make block diagonal pi;
  bPi <- blockLambda_fcn(pi, p, D)

  
  theta <- solve(t(kron) %*% bPi %*%
                       scaled_XTX %*% kron) %*%
    t(kron) %*% bPi %*%
    scaled_XTX %*% beta

  theta
}



## Asymptotic confidence intervals:
kld_linear_aCI <- function(pi, scaled_XTX, indv_mle, p, D, kron, alpha){
  # Make block diagonal pi;
  if (all.equal(pi, rep(0, D))==T){var <- solve(scaled_XTX)
  } else{
  bPi <- blockLambda_fcn(pi, p, D)
  I <- diag(1, nrow = p*D)
  A <- solve(t(kron) %*% bPi %*% scaled_XTX %*% kron) %*% t(kron) %*% bPi %*% scaled_XTX

  
  
  theta_hat <- A %*% indv_mle

  beta_hat <- (I-bPi) %*% indv_mle +
    bPi %*% kron %*%
    theta_hat

  beta_hat_prime <- I - bPi +  bPi %*% kron %*% A
  interior <- solve(scaled_XTX)
  var <- beta_hat_prime %*% interior %*% t(beta_hat_prime)
  }

  z_crit <- qnorm(1-alpha/2, mean = 0, sd = 1)
  lb <- beta_hat - z_crit * sqrt(diag(var))
  ub <- beta_hat + z_crit * sqrt(diag(var))

  list(lb = lb, ub = ub)
}



#### bbeta_m ci ######
ham_bbetam_aci <- function(pi, scaled_XTX, tau, indv_mle, p, D, kron, alpha){
  # Make block diagonal pi;
  bPi <- blockLambda_fcn(pi, p, D)

  
  
  theta_hat <- solve(t(kron) %*% bPi %*%
                       scaled_XTX %*% kron) %*%
    t(kron) %*% bPi %*% scaled_XTX %*% indv_mle

  beta_hat <- (diag(1, nrow = p*D)-bPi) %*% indv_mle +
    bPi %*% kron %*%
    theta_hat

  hbbeta_m <- 1/D*t(kron(beeta_hat))

  beta_hat_prime <- diag(1, nrow = p*D) - bPi +
    bPi %*% kron %*% solve(t(kron) %*% bPi %*%
                             scaled_XTX %*% kron) %*%
    t(kron) %*% bPi %*% scaled_XTX
  interior <- solve(scaled_XTX)
  var <- beta_hat_prime %*% interior %*% t(beta_hat_prime)

  z_crit <- qnorm(1-alpha/2, mean = 0, sd = 1)
  lb <- beta_hat - z_crit * sqrt(diag(var))
  ub <- beta_hat + z_crit * sqrt(diag(var))

  list(lb = lb, ub = ub)

}



