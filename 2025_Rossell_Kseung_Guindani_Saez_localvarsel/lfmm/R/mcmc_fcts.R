


#' The Half Cauchy Distribution
#' 
#' Probability density function for the half Cauchy distribution with scale 
#' equal to sigma.
#' 
#' @param x Vector of scalars.
#' @param sigma Scale parameter.
#' @return The pdf of the half Cauchy distribution evaluated at points x.
dhalfcauhy <- function(x, sigma, log = T){
  
  out <- log(2) + log(sigma) - log(pi) - log(sigma^2 + x^2)
  if (!log){
    out <- exp(out)
  }
  
  return(out)
}


#' Smoothness Inducing Penalty
#' 
#' Construct the covariance matrix P of the smoothness inducing prior for the
#' spline coefficients.
#' 
#' @param K Number of knots.
#' @return Covariance of the smoothness inducing prior (penalizing first 
#'         differences in beta).
P_smooth1 <- function(K){
  
  D <- diag(rep(1,K))
  D <- diff(D)
  P <- crossprod(D)
  
  return(P)
}


#' B-spline Basis
#' 
#' Generate the B-spline basis matrix using De Boor's recursive formula.
#' 
#' @param xgrid Vector of scalars where to evaluate the spline basis.
#' @param knots Vector of knot points.
#' @param q Spline basis degree.
#' @return A matrix of dimension (length(xgrid) \times length(knots) + q - 1).
B_basis_q <- function(xgrid, knots, q = 1){
  n <- length(xgrid) # number of grid points where to evaluate the spline
  K <- length(knots) # number of knots
  
  B_k <- list()
  
  B_k[[1]] <- array(0, dim = c(n, K - 1))
  for (j in 1:(K - 1)){
    if (j < (K - 1)){
      B_k[[1]][which((xgrid >= knots[j]) & (xgrid < knots[j+1])),j] <- 1
    }
    else if (j == (K - 1)){
      B_k[[1]][which((xgrid >= knots[j]) & (xgrid <= knots[j+1])),j] <- 1
    }
  }
  B_temp <- cbind(rep(0, n), B_k[[1]], rep(0, n))
  knots <- c(knots[1], knots, knots[length(knots)])
  
  for (k in 1:q){
    B_k[[k+1]] <- array(0, dim = c(n, K - 1 + k))
    for (j in 1:(K -1 + k)){
      if (j == 1){
        B_k[[k+1]][,j] <- (knots[j+k+1] - xgrid)/(knots[j+k+1] - knots[j+1]) * B_temp[,j+1]
        
      }
      else if (j == (K -1 + k)){
        B_k[[k+1]][,j] <- (xgrid - knots[j])/(knots[j+k] - knots[j]) * B_temp[,j]
        
      }
      else {
        B_k[[k+1]][,j] <- (xgrid - knots[j])/(knots[j+k] - knots[j]) * B_temp[,j] + (knots[j+k+1] - xgrid)/(knots[j+k+1] - knots[j+1]) * B_temp[,j+1]
        
      }
    }
    B_temp <- cbind(rep(0, n), B_k[[k+1]], rep(0, n))
    knots <- c(knots[1], knots, knots[length(knots)])
  }
  
  return(B_k[[q + 1]])
}


#' Hamming Ball
#' 
#' Computes the Hamming Ball centered at x with radius r.
#' 
#' @param x Center of the Hamming Ball.
#' @param S Number of states .
#' @param r Radius of the Hamming Ball.
#' @return Hamming Ball
H_ball <- function(x, S, r){
  
  K <- length(x) # number of chains
  card_m <- (S - 1)^(0:r) * choose(K, 0:r) 
  # sum(card_m) is the cardinality of each HB with radius r
  
  # First, put the current vector in its HB
  HB_temp <- matrix(x, K, 1)
  
  for (m in 1:r){ # loop the possible HB radius
    HB_diff <- matrix(rep(x, card_m[m+1]), K, card_m[m+1])
    index <- utils::combn(K, m) # what elements of the vector to change?
    
    for (j in 1:ncol(index)){
      vec <- NULL
      for (i in 1:nrow(index)){ 
        vec <- c(vec, (1:S)[-x[index[i,j]]])
      }
      prop_col <- t(gtools::permutations(n = length(unique(vec)), r = m, v = vec, repeats.allowed = TRUE))
      keep_col <- prop_col[,which(colSums(prop_col != x[index[,j]]) == m)]
      HB_diff[index[,j], ((j-1)*(card_m[m+1]/ncol(index))+1):(j*card_m[m+1]/ncol(index))] <- keep_col
    }
    HB_temp <- cbind(HB_temp,HB_diff) # save these in the Hamming Ball
  }
  
  return (HB_temp)
}


#' Hamming Ball
#' 
#' Samples a configuration uniformly within the HB
#' 
#' @param x Center of the Hamming Ball.
#' @param S Number of states .
#' @param r Radius of the Hamming Ball.
#' @return uniform sample within the Hamming Ball
H_ball_unif <- function(x, S, r){
  HB_temp <- H_ball(x, S, r)
  HB_samp <- HB_temp[,sample(1:ncol(HB_temp), 1)]
  
  return (HB_samp)
}



#' Local Clustering with Random Effects
#' 
#' Markov Chain Monte Carlo for the local clustering model with a simple 
#' categorical predictor (case with random effects). 
#' 
#' @param y Vector of observations.
#' @param time Vector of time points corresponding to the observations.
#' @param ind Vector of participants indicators corresponding to the observations.
#' @param X Matrix of categorical predictors.
#' @param xgrid Vector where to evaluate the posterior draws for the functional 
#'              parameters.
#' @param hypers List containing the hyperparameters $a_{\sigma}, b_{\sigma}, 
#'               a_{\alpha}, b_{\alpha}, \sigma_{s}$.
#' @param Niter Number of iterations for the MCMC. 
#' @param burnin Burnin for the MCMC. 
#' @param thin Thinning factor for the MCMC. 
#' @return A list of MCMC objects.
loc_clust_reff <- function(y, time, ind, X, xgrid, hypers, Niter = 6000, 
                           burnin = 1000, thin = 1){
  
  # Auxiliary variables
  # std <- sd(y)
  std <- 1
  y <- y / std
  samp_size <- (Niter - burnin)/thin # sample size 
  n <- length(y) # number of observations
  n_ind <- length(unique(ind))
  n_grid <- length(xgrid)
  T_min <- min(time)
  T_max <- max(time)
  knots <- T_min:T_max
  time <- time - T_min + 1
  K <- length(knots)
  knots <- 1:K
  xgrid <- xgrid - T_min + 1
  d_1 <- length(unique(X)) # number of levels for x
  B <- B_basis_q(time, knots)
  Bgrid <- B_basis_q(xgrid, knots)
  P_u <- P_smooth1(K)
  
  # Initialize mu_0, sigma_0
  mu_0 <- aggregate(y, list(time), mean)$x
  sigma2_0 <- aggregate(y, list(time), var)$x
  
  # Set MCMC objects
  post_mean <- array(NA, dim = c(n_grid, d_1, samp_size))
  z_chain <- array(NA, dim = c(d_1, K, samp_size))
  post_data <- array(NA, dim = c(samp_size, n))
  sigma2_y_chain <- array(NA, dim = samp_size)
  sigma2_1_chain <- array(NA, dim = samp_size)
  alpha_chain <- array(NA, dim = samp_size)
  post_reff_gr <- array(NA, dim = c(n_grid, n_ind, samp_size))
  sigma2_us_chain <- array(NA, dim = samp_size)
  sigma2_ua_chain <- array(NA, dim = samp_size)
  
  # Set initial values
  sigma2_y_old <- 10
  sigma2_1_old <- 1
  nu_old <- 1
  alpha_old <- 1
  z_old <- array(NA, dim = c(d_1, K))
  for (h in 1:d_1){
    z_old[h,] <- h
  }
  B_beta <- B_beta_u <- array(0, n)
  beta_old <- beta_star_old <- array(NA, dim = c(d_1, K))
  for (h in 1:d_1){
    for (k in 1:K){
      beta_old[h,k] <- beta_star_old[h,k] <- mean(y[which((time == k) & (X == h))])
    }
    B_beta[X == h] <- B[X == h,] %*% beta_old[h,] 
  }
  pi_0_old <- rdirichlet(alpha_old/d_1 + table_int(z_old[,1], d_1))
  Q_old <- array(NA, dim = c(d_1, d_1))
  tr_count <- count_assign(z_old)
  for (h in 1:d_1){
    Q_old[h,] <- rdirichlet(alpha_old/d_1 + tr_count[h,])
  }
  sigma2_us_old <- 0.5
  sigma2_ua_old <- 0.5
  
  # Auxiliary quantities
  beta_mess <- array(NA, dim = c(d_1, K))
  sd_MH_alpha <- sd_MH_sigma <- sd_MH_sigma_ua <- sd_MH_sigma_us <- 0.1
  acc_alpha <- acc_sigma <- acc_sigma_ua <- acc_sigma_us <- 0
  n_batch <- 0
  
  # Gibbs Sampler
  it <- 1
  pb <- txtProgressBar(style = 3)
  for (iter in 1:Niter){
    # Adaptive MH for alpha
    if (iter %% 20 == 0){
      n_batch <- n_batch + 1
      delta_n <- min(0.01, n_batch^(-0.5))
      
      if (acc_alpha/iter > 0.44){
        sd_MH_alpha <- exp(log(sd_MH_alpha) + delta_n)
      }
      else{
        sd_MH_alpha <- exp(log(sd_MH_alpha) - delta_n)
      }
      if (acc_sigma/iter > 0.44){
        sd_MH_sigma <- exp(log(sd_MH_sigma) + delta_n)
      }
      else{
        sd_MH_sigma <- exp(log(sd_MH_sigma) - delta_n)
      }
      if (acc_sigma_ua/iter > 0.44){
        sd_MH_sigma_ua <- exp(log(sd_MH_sigma_ua) + delta_n)
      }
      else{
        sd_MH_sigma_ua <- exp(log(sd_MH_sigma_ua) - delta_n)
      }
      if (acc_sigma_us/iter > 0.44){
        sd_MH_sigma_us <- exp(log(sd_MH_sigma_us) + delta_n)
      }
      else{
        sd_MH_sigma_us <- exp(log(sd_MH_sigma_us) - delta_n)
      }
    }
    
    # Compute residuals: data minus random effects
    res_y <- y - B_beta_u
    
    # (1) Update the cluster specific curves (dynamic updating)
    for (k in 1:K){ # loop over spline coefficients
      for (h in 1:d_1){ # indicator for the latent
        X_1k <- which(z_old[,k] == h) # X_{1,k}
        if (length(X_1k) > 0){ # h \in \mathcal{Z}_{1,k}: posterior update
          
          # Pick data with covariate levels clustered in group h and at the
          # correct time stamps
          idx_i <- which( (X %in% X_1k) & (time == k) )
          n_all <- length(idx_i)
          y_temp <- res_y[idx_i]
          X_i <- X[which((X %in% X_1k) & (time == k))]
          
          
          var_dat <- sigma2_y_old/n_all
          m_dat <- sum(y_temp)/n_all
          if (k == 1){
            var_smooth <- sigma2_1_old/length(unique(z_old[X_1k,k+1]))
          }
          else if (k == K){
            var_smooth <- sigma2_1_old/length(unique(z_old[X_1k,k-1]))
          }
          else {
            var_smooth <- sigma2_1_old/length(c(unique(z_old[X_1k,k-1]), unique(z_old[X_1k,k+1])))
          }
          if (k == 1){
            m_smooth <- sum(beta_star_old[unique(z_old[X_1k,k+1]),k+1])/length(unique(z_old[X_1k,k+1]))
          }
          else if (k == K){
            m_smooth <- sum(beta_star_old[unique(z_old[X_1k,k-1]),k-1])/length(unique(z_old[X_1k,k-1]))
          }
          else {
            m_smooth <- (sum(beta_star_old[unique(z_old[X_1k,k-1]),k-1]) + sum(beta_star_old[unique(z_old[X_1k,k+1]),k+1]))/length(c(unique(z_old[X_1k,k-1]), unique(z_old[X_1k,k+1])))
          }
          var_post <- 1/(1/var_dat + 1/var_smooth)
          m_post <- var_post * (m_dat/var_dat + m_smooth/var_smooth)
          
          beta_star_old[h,k] <- rnorm(1, m_post, sqrt(var_post))
        }
        else { # h \notin \mathcal{Z}_{1,k}: prior sampling
          beta_star_old[h,k] <- rnorm(1, mu_0, sqrt(sigma2_0))
        }
      }
    }
    
    
    # (2) Update the cluster assignments
    for (x1 in 1:d_1){
      
      # (a) Pass messages backwards
      
      # Final condition
      beta_mess[,K] <- 1/d_1
      for (k in (K - 1):1){
        idx <- which( (X == x1) & (time == (k+1)))
        
        prob_mat <- backprop_mess_lin(res_y[idx], beta_star_old[,k:(k+1)], 
                                      sigma2_y_old, beta_mess[,k+1], Q_old, 
                                      sigma2_1_old)
        prob_vec <- sum_rows_log(prob_mat)
        
        # Normalize messages for stability
        beta_mess[,k] <- as.numeric(normalize_log(prob_vec))
      }
      
      # (b) Sample states forward
      
      # Initial condition
      idx <- which( (X == x1) & (time == 1))
      prob_vec <- as.numeric(log(pi_0_old) + 
                               log(beta_mess[,1]) - 
                               0.5/sigma2_y_old * colSums((outer(res_y[idx], beta_star_old[,1], '-'))^2))
      
      # Normalize probabilities
      prob_vec <- as.numeric(normalize_log(prob_vec))
      z_old[x1,1] <- sample(1:d_1, 1, FALSE, prob_vec)
      
      for (k in 2:K){
        idx <- which( (X == x1) & (time == k))
        prob_vec <- as.numeric(log(Q_old[z_old[x1,k-1],]) + 
                                 log(beta_mess[,k]) - 
                                 0.5 / sigma2_1_old * ((beta_star_old[,k] - beta_star_old[z_old[x1,k-1],k-1])^2) - 
                                 0.5/sigma2_y_old * colSums((outer(res_y[idx], beta_star_old[,k], '-'))^2))
        
        
        # Normalize probabilities
        prob_vec <- as.numeric(normalize_log(prob_vec))
        z_old[x1,k] <- sample(1:d_1, 1, FALSE, prob_vec)
      }
    }
    
    
    # (3) Update the transition probabilities
    # Use a different dynamic for each covariate
    pi_0_old <- rdirichlet(alpha_old/d_1 + table_int(z_old[,1], d_1))
    tr_count <- count_assign(z_old)
    for (h in 1:d_1){
      Q_old[h,] <- rdirichlet(alpha_old/d_1 + tr_count[h,])
    }
    # Prevent from underflow
    pi_0_old[pi_0_old == 0] <- 1.0e-323
    Q_old[Q_old == 0] <- 1.0e-323
    
    
    # (4) Update the hyperparameter for the transition probabilities
    alpha_temp <- exp(rnorm(1, log(alpha_old), sd_MH_alpha))
    acc_rate <- min(c(0, (d_1 + 1) * (lgamma(alpha_temp) - d_1 * lgamma(alpha_temp/d_1)) + 
                        (alpha_temp/d_1 - 1) * (sum(log(Q_old)) + sum(log(pi_0_old))) +
                        dgamma(alpha_temp, hypers$a_alpha, hypers$b_alpha, log = T) +
                        log(alpha_temp) -
                        (d_1 + 1) * (lgamma(alpha_old) - d_1 * lgamma(alpha_old/d_1)) - 
                        (alpha_old/d_1 - 1) * (sum(log(Q_old)) + sum(log(pi_0_old))) -
                        dgamma(alpha_old, hypers$a_alpha, hypers$b_alpha, log = T) -
                        log(alpha_old)))
    lu <- log(runif(1))
    if (lu < acc_rate){
      alpha_old <- alpha_temp
      acc_alpha <- acc_alpha + 1
    }
    
    
    # (5) Assign the cluster specific curves
    for (h in 1:d_1){
      beta_old[h,] <- beta_star_old[cbind(z_old[h,],1:K)]
      B_beta[X == h] <- B[X == h,] %*% beta_old[h,]
    }
    
    
    # (6) Update the cluster specific smoothness parameters
    RSS <- 0
    for (h in 1:d_1){
      RSS <- RSS + as.numeric(crossprod(beta_old[h,], P_u) %*% beta_old[h,])
    }
    sigma2_1_old <- 1/rgamma(1, 0.5 + 0.5 * (d_1*K),
                             1/nu_old + 0.5 * RSS)
    nu_old <- 1/rgamma(1, 0.5, 1/sigma2_1_old + 1/hypers$sigma_s^2)
    # sigma2_1_temp <- exp(rnorm(1, log(sigma2_1_old), sd_MH_sigma))
    # acc_rate <- min(c(0, -0.5 * (K * d_1) * log(sigma2_1_temp) - 
    #                     0.5/sigma2_1_temp * RSS + 
    #                     dhalfcauhy(sqrt(sigma2_1_temp), hypers$sigma_s, log = T) + 
    #                     log(sigma2_1_temp) +
    #                     0.5 * (K * d_1) * log(sigma2_1_old) + 
    #                     0.5/sigma2_1_old * RSS -
    #                     dhalfcauhy(sqrt(sigma2_1_old), hypers$sigma_s, log = T) -
    #                     log(sigma2_1_old)))
    # lu <- log(runif(1))
    # if (lu < acc_rate){
    #   sigma2_1_old <- sigma2_1_temp
    #   acc_sigma <- acc_sigma + 1
    # }
    
    
    # Compute residuals: data minus main effects
    res_u <- y - B_beta
    
    # (7) Update random effects
    reff <- sample_reff(res_u, B, P_u, ind, sigma2_y_old, sigma2_us_old, 
                        sigma2_ua_old)
    beta_ui_old <- reff$beta_ui
    B_beta_u <- as.numeric(reff$B_beta_u)
    
    
    # (8) Update random effects variances
    sigma2_ua_prop <- exp(rnorm(1, mean = log(sigma2_ua_old), sd = sd_MH_sigma_ua))
    acc_rate <- min(c(0, 0.5 * n_ind * log(det(diag(K)/sigma2_ua_prop + P_u/sigma2_us_old)) -
                        0.5 * sum(diag(beta_ui_old %*% tcrossprod(diag(K)/sigma2_ua_prop, beta_ui_old))) +
                        dhalfcauhy(sqrt(sigma2_ua_prop), hypers$sigma_s, TRUE) + 
                        log(sigma2_ua_prop) -
                        0.5 * n_ind * log(det(diag(K)/sigma2_ua_old + P_u/sigma2_us_old)) +
                        0.5 * sum(diag(beta_ui_old %*% tcrossprod(diag(K)/sigma2_ua_old, beta_ui_old))) -
                        dhalfcauhy(sqrt(sigma2_ua_old), hypers$sigma_s, TRUE) - 
                        log(sigma2_ua_old)))
    lu <- log(runif(1))
    if (lu < acc_rate){
      sigma2_ua_old <- sigma2_ua_prop
      acc_sigma_ua <- acc_sigma_ua + 1
    }
    sigma2_us_prop <- exp(rnorm(1, mean = log(sigma2_us_old), sd = sd_MH_sigma_us))
    acc_rate <- min(c(0, 0.5 * n_ind * log(det(diag(K)/sigma2_ua_old + P_u/sigma2_us_prop)) -
                        0.5 * sum(diag(beta_ui_old %*% tcrossprod(P_u/sigma2_us_prop, beta_ui_old))) +
                        dhalfcauhy(sqrt(sigma2_us_prop), hypers$sigma_s, TRUE) + 
                        log(sigma2_us_prop) -
                        0.5 * n_ind * log(det(diag(K)/sigma2_ua_old + P_u/sigma2_us_old)) +
                        0.5 * sum(diag(beta_ui_old %*% tcrossprod(P_u/sigma2_us_old, beta_ui_old))) -
                        dhalfcauhy(sqrt(sigma2_us_old), hypers$sigma_s, TRUE) - 
                        log(sigma2_us_old)))
    if (lu < acc_rate){
      sigma2_us_old <- sigma2_us_prop
      acc_sigma_us <- acc_sigma_us + 1
    }
    
    
    # (9) Update the variance of the data (global)
    RSS <- as.numeric(crossprod(y - B_beta - B_beta_u))
    sigma2_y_old <- 1/rgamma(1, hypers$a_sigma + n/2, hypers$b_sigma + 0.5 * RSS)
    
    
    # Save parameters in the chain
    if ( (iter > burnin) && (iter %% thin == 0) ){
      for (h in 1:d_1){
        post_mean[,h,it] <- std * Bgrid %*% beta_old[h,]
        z_chain[h,,it] <- z_old[h,]
      }
      for (i in 1:n_ind){
        post_reff_gr[,i,it] <- std * Bgrid %*% beta_ui_old[i,]
      }
      post_data[it,] <- std * (B_beta + B_beta_u)
      sigma2_1_chain[it] <- sigma2_1_old
      sigma2_y_chain[it] <- sigma2_y_old
      alpha_chain[it] <- alpha_old
      sigma2_us_chain[it] <- sigma2_us_old
      sigma2_ua_chain[it] <- sigma2_ua_old
      
      it <- it + 1
    }
    setTxtProgressBar(pb, iter/Niter)
  }
  
  return(list('Z' = z_chain, 'post_mean' = post_mean, 'post_data' = post_data,
              'post_reff_gr' = post_reff_gr, 'sigma2_us' = sigma2_us_chain, 
              'sigma2_ua' = sigma2_ua_chain, 'sigma2_1' = sigma2_1_chain, 
              'sigma2_y' = sigma2_y_chain, 'alpha' = alpha_chain))
}



#' Local Clustering of multiple predictors with Random Effects
#' 
#' Markov Chain Monte Carlo for the local clustering model with multiple 
#' categorical predictors (case with random effects). 
#' 
#' @param y Vector of observations.
#' @param time Vector of time points corresponding to the observations.
#' @param ind Vector of participants indicators corresponding to the observations.
#' @param X Matrix of categorical predictors.
#' @param xgrid Vector where to evaluate the posterior draws for the functional 
#'              parameters.
#' @param hypers List containing the hyperparameters $a_{\sigma}, b_{\sigma}, 
#'               a_{\alpha}, b_{\alpha}, \sigma_{s}$.
#' @param Niter Number of iterations for the MCMC. 
#' @param burnin Burnin for the MCMC. 
#' @param thin Thinning factor for the MCMC. 
#' @return A list of MCMC objects.
loc_clust_multi_reff <- function(y, time, ind, X, t_grid, Xgrid, hypers, 
                                 Niter = 6000, burnin = 1000, thin = 1){
  
  # Auxiliary variables
  n_ind <- length(unique(ind))
  samp_size <- (Niter - burnin)/thin # sample size 
  n <- length(y) # number of observations
  T_min <- min(time)
  T_max <- max(time)
  knots <- T_min:T_max
  time <- time - T_min + 1
  K <- length(knots)
  p <- ncol(X)
  d_j <- apply(X, 2, function(x){length(unique(x))}) # number of levels for x
  B <- B_basis_q(time, knots)
  Bgrid <- B_basis_q(t_grid, knots)
  P_u <- P_smooth1(K)
  M_max <- min(5, floor(sqrt(prod(d_j))))
  max_cov <- 5
  
  # Initialize mu_0, sigma_0
  mu_0 <- aggregate(y, list(time), mean)$x
  sigma2_0 <- aggregate(y, list(time), var)$x
  
  # Set MCMC objects
  phi_chain <- array(NA, dim = c(samp_size, p))
  l_k_chain <- array(NA, dim = c(samp_size, p, K))
  post_pred <- array(NA, dim = c(samp_size, nrow(Xgrid), K))
  post_data <- array(NA, dim = c(samp_size, n))
  post_mean <- array(NA, dim = c(samp_size, n))
  sigma2_y_chain <- array(NA, dim = samp_size)
  sigma2_1_chain <- array(NA, dim = samp_size)
  alpha_chain <- array(NA, dim = c(samp_size, p))
  RSS_chain <- array(NA, dim = samp_size)
  sigma2_us_chain <- array(NA, dim = samp_size)
  sigma2_ua_chain <- array(NA, dim = samp_size)
  post_reff_gr <- array(NA, dim = c(length(t_grid), n_ind, samp_size))
  
  # Set initial values
  phi_old <- rep(1, p)
  l_k_old <- array(1, dim = c(p, K))
  l_k_prop <- array(1, dim = c(p, K))
  sigma2_y_old <- 10
  sigma2_1_old <- 1
  sigma2_us_old <- 0.1
  sigma2_ua_old <- 0.1
  nu_old <- 1
  alpha_old <- rep(1, p)
  alpha_2_old <- 1
  z_old <- list()
  for (j in 1:p){
    z_old[[j]] <- array(1, dim = c(d_j[j], K))
  }
  z_sl_old <- list()
  for (k in 1:K){
    joint_part <- as.matrix(expand.grid(lapply(z_old, function(x){unique(x[,k])})))
    z_sl_old[[k]] <- 1:nrow(joint_part)
  }
  
  
  B_beta <- array(0, n)
  B_beta_pred <- array(0, dim = c(nrow(Xgrid), K))
  B_beta_u <- rep(0, n)
  beta_star_old <- list()
  for (k in 1:K){
    beta_star_old[[k]] <- array(NA, dim = l_k_old[,k])
    cl_cov <- expand.grid(lapply(l_k_old[,k], function(x){1:x}))
    for (h in 1:nrow(cl_cov)){
      idx_cov <- rep(T, n)
      for (j in 1:p){
        idx_cov <- idx_cov & (X[,j] %in% which(z_old[[j]][,k] == cl_cov[h,j]))
      }
      idx_i <- which((time == k) & idx_cov)
      beta_star_old[[k]][matrix(as.numeric(cl_cov[h,]), 1, p)] <- mean(y[idx_i])
      B_beta[idx_i] <- beta_star_old[[k]][matrix(as.numeric(cl_cov[h,]), 1, p)]
    }
  }
  
  pi_0_old <- Q_old <- list()
  for (j in 1:p){
    pi_0_old[[j]] <- rdirichlet(alpha_old[j]/d_j[j] + table_int(z_old[[j]][,1], d_j[j]))
    Q_old[[j]] <- array(NA, dim = c(d_j[j], d_j[j]))
    tr_count <- count_assign(z_old[[j]])
    for (h in 1:d_j[j]){
      Q_old[[j]][h,] <- rdirichlet(alpha_old[j]/d_j[j] + tr_count[h,])
    }
  }
  pi_2_old <- rep(1/M_max, M_max)
  
  # Auxiliary quantities
  sd_MH_alpha <- rep(0.05, p)
  sd_MH_phi <- rep(0.05, p)
  sd_MH_sigma <- 0.1
  sd_MH_sigma_ua <- 0.1
  sd_MH_sigma_us <- 0.1
  acc_alpha <- rep(0, p)
  acc_phi <- rep(0, p)
  acc_sigma <- 0
  acc_sigma_ua <- 0
  acc_sigma_us <- 0
  n_batch <- 0
  
  # Gibbs Sampler
  it <- 1
  pb <- txtProgressBar(style = 3)
  for (iter in 1:Niter){
    # Adaptive MH for alpha
    if (iter %% 50 == 0){
      n_batch <- n_batch + 1
      delta_n <- min(0.01, n_batch^(-0.5))
      
      for (j in 1:p){
        if (acc_alpha[j]/iter > 0.44){
          sd_MH_alpha[j] <- exp(log(sd_MH_alpha[j]) + delta_n)
        }
        else{
          sd_MH_alpha[j] <- exp(log(sd_MH_alpha[j]) - delta_n)
        }
        if (acc_phi[j]/iter > 0.44){
          sd_MH_phi[j] <- exp(log(sd_MH_phi[j]) + delta_n)
        }
        else{
          sd_MH_phi[j] <- exp(log(sd_MH_phi[j]) - delta_n)
        }
      }
      if (acc_sigma/iter > 0.44){
        sd_MH_sigma <- exp(log(sd_MH_sigma) + delta_n)
      }
      else{
        sd_MH_sigma <- exp(log(sd_MH_sigma) - delta_n)
      }
      if (acc_sigma_ua/iter > 0.44){	
        sd_MH_sigma_ua <- exp(log(sd_MH_sigma_ua) + delta_n)	
      }	
      else{	
        sd_MH_sigma_ua <- exp(log(sd_MH_sigma_ua) - delta_n)	
      }	
      if (acc_sigma_us/iter > 0.44){	
        sd_MH_sigma_us <- exp(log(sd_MH_sigma_us) + delta_n)	
      }	
      else{	
        sd_MH_sigma_us <- exp(log(sd_MH_sigma_us) - delta_n)	
      }
    }
    
    # Compute residuals: data minus random effects	
    res_y <- y - B_beta_u
    
    RSS <- 0
    n_RSS <- 0
    for (k in 1:K){ # loop over spline coefficients
      
      # (1) Update the cluster labels
      y_temp <- res_y[which(time == k)]
      X_temp <- X[which(time == k),]
      l_k_prop[,k] <- l_k_old[,k]
      
      for (j in 1:p){ # loop over covariates
        
        if (sum(l_k_old[,k] > 1) <= max_cov){
          
          # Case for generic r = 1 (slower implementation)
          # z_prop <- H_ball_unif(z_old[[j]][,k], d_j[j], 1)
          
          # joint_part
          # z_sl_old[[k]]
          
          # Special case when r = 1 (faster implementation)
          z_prop <- z_old[[j]][,k]
          z_prop[sample(1:d_j[j], 1)] <- sample(1:d_j[j], 1)
          
          if ( sum(z_prop != z_old[[j]][,k]) > 0 ){
            
            l_k_prop[j,k] <- length(unique(z_prop))
            
            z_prop_list <- z_old
            z_prop_list[[j]][,k] <- z_prop
            
            if (k == 1){
              pr_old <- sum(log(pi_0_old[[j]][z_old[[j]][,k]])) +
                sum(log(Q_old[[j]][cbind(z_old[[j]][,k], z_old[[j]][,k+1])]))
              
              pr_prop <- sum(log(pi_0_old[[j]][z_prop])) +
                sum(log(Q_old[[j]][cbind(z_prop, z_old[[j]][,k+1])]))
            }
            else if (k == K){
              pr_old <- sum(log(Q_old[[j]][cbind(z_old[[j]][,k-1], z_old[[j]][,k])]))
              
              pr_prop <- sum(log(Q_old[[j]][cbind(z_old[[j]][,k-1], z_prop)]))
            }
            else {
              pr_old <- sum(log(Q_old[[j]][cbind(z_old[[j]][,k-1], z_old[[j]][,k])])) +
                sum(log(Q_old[[j]][cbind(z_old[[j]][,k], z_old[[j]][,k+1])]))
              
              pr_prop <- sum(log(Q_old[[j]][cbind(z_old[[j]][,k-1], z_prop)])) +
                sum(log(Q_old[[j]][cbind(z_prop, z_old[[j]][,k+1])]))
            }
            
            # Explicitly compute the implied joint partition by the
            # marginals
            joint_part <- as.matrix(expand.grid(lapply(z_old, function(x){unique(x[,k])})))
            joint_part_prop <- as.matrix(expand.grid(lapply(z_prop_list, function(x){unique(x[,k])})))
            z_sl_prop <- sample(1:M_max, nrow(joint_part_prop), replace = T, prob = pi_2_old)
            
            lu <- log(runif(1))
            acc_rate <- min(c(0, 
                              log_lik_z(y_temp, X_temp, z_prop_list, sigma2_y_old, mu_0[k], sigma2_0[k], joint_part_prop, z_sl_prop, k) +
                                pr_prop + 
                                sum(log(pi_2_old[z_sl_prop])) -
                                phi_old[j] * l_k_prop[j,k] - (
                                  log_lik_z(y_temp, X_temp, z_old, sigma2_y_old, mu_0[k], sigma2_0[k], joint_part, z_sl_old[[k]], k) +
                                    pr_old + 
                                    sum(log(pi_2_old[z_sl_old[[k]]])) -
                                    phi_old[j] * l_k_old[j,k])))
            if (lu < acc_rate){
              z_old <- z_prop_list
              l_k_old[j,k] <- l_k_prop[j,k]
              z_sl_old[[k]] <- z_sl_prop
            }
          }
        }
      }
      
      
      # (2) Update the cluster specific curves (dynamic updating)
      beta_star_old[[k]] <- array(NA, dim = l_k_old[,k])
      joint_part <- as.matrix(expand.grid(lapply(z_old, function(x){unique(x[,k])})))
      cl_cov <- expand.grid(lapply(l_k_old[,k], function(x){1:x}))
      
      
      z_sl_un <- unique(z_sl_old[[k]])
      for (h in z_sl_un){ # indicator for the latent
        idx_z <- which(z_sl_old[[k]] == h)
        idx_cov <- rep(0, n)
        idx_cov_pred <- rep(0, nrow(Xgrid))
        cov_hh <- array(NA, dim = c(length(idx_z), p))
        for (hh in 1:length(idx_z)){
          idx_cov <- idx_cov + find_cov(X, z_old, as.numeric(joint_part[idx_z[hh],]), k)
          idx_cov_pred <- idx_cov_pred + find_cov(Xgrid, z_old, as.numeric(joint_part[idx_z[hh],]), k)
          # cov_hh[hh,] <- find_cov_2(z_old, as.numeric(joint_part[idx_z[hh],]), k)
        }
        # as.numeric(beta_star_old[[k-1]])
        # as.numeric(beta_star_old[[k+1]])
        
        
        # Pick data with covariate levels clustered in group h and at the
        # correct time stamps
        idx_i <- which((time == k) & idx_cov)
        idx_i_pred <- which(as.numeric(idx_cov_pred) == 1)
        
        n_all <- length(idx_i)
        y_temp <- res_y[idx_i]
        X_i <- X[idx_i,]

        var_dat <- sigma2_y_old/n_all
        m_dat <- sum(y_temp)/n_all          
        
        # if (k == 1){
        #   m_post
        #   n_post
        #   var_smooth <- sigma2_1_old / n_post
        #   m_smooth <- sum(m_post) / n_post
        #   # var_smooth <- sigma2_1_old/length(unique(z_old[X_1k,k+1]))
        #   # m_smooth <- sum(beta_star_old[[k+1]][unique(z_old[X_1k,k+1]),k+1])/length(unique(z_old[X_1k,k+1]))
        # }
        # else if (k == K){
        #   m_pre
        #   n_pre
        #   var_smooth <- sigma2_1_old / n_pre
        #   m_smooth <- sum(m_pre) / n_pre
        #   # var_smooth <- sigma2_1_old/length(unique(z_old[X_1k,k-1]))
        #   # m_smooth <- sum(beta_star_old[[k-1]][unique(z_old[X_1k,k-1]),k-1])/length(unique(z_old[X_1k,k-1]))
        # }
        # else {
        #   m_pre
        #   m_post
        #   n_pre
        #   n_post
        #   var_smooth <- sigma2_1_old / (n_pre + n_post)
        #   m_smooth <- sum(c(m_pre, m_post)) / (n_pre + n_post)
        #   
        #   # var_smooth <- sigma2_1_old/length(c(unique(z_old[X_1k,k-1]), unique(z_old[X_1k,k+1])))
        #   # m_smooth <- (sum(beta_star_old[[k-1]][unique(z_old[X_1k,k-1]),k-1]) + sum(beta_star_old[[k+1]][unique(z_old[X_1k,k+1]),k+1]))/length(c(unique(z_old[X_1k,k-1]), unique(z_old[X_1k,k+1])))
        # }
        
        ## REMOVE
        m_smooth <- 0
        var_smooth <- +Inf
        ## REMOVE 

        var_star <- 1/(1/var_dat + 1/var_smooth)
        m_star <- var_star * (m_dat/var_dat + m_smooth/var_smooth)
        
        beta_temp <- rnorm(1, m_star, sqrt(var_star))
        # if (k > 1){
        #   RSS <- RSS + (beta_temp - m_pre)^2
        #   n_RSS <- n_RSS + n_pre
        # }
        
        for (hh in 1:length(idx_z)){
          beta_star_old[[k]][matrix(as.numeric(cl_cov[idx_z[hh],]), 1, p)] <- beta_temp
        }
        
        B_beta[idx_i] <- beta_star_old[[k]][matrix(as.numeric(cl_cov[idx_z[1],]), 1, p)]
        B_beta_pred[idx_i_pred,k] <- beta_star_old[[k]][matrix(as.numeric(cl_cov[idx_z[1],]), 1, p)]
      }
    }
    
    
    # (3) Update the transition probabilities
    # Use a different dynamic for each covariate
    for (j in 1:p){
      pi_0_old[[j]] <- rdirichlet(alpha_old[j]/d_j[j] + table_int(z_old[[j]][,1], d_j[j]))
      for (h in 1:d_j[j]){
        tr_count <- count_assign(z_old[[j]])
        Q_old[[j]][h,] <- rdirichlet(alpha_old[j]/d_j[j] + tr_count[h,])
      }
      # Prevent from underflow
      pi_0_old[[j]][pi_0_old[[j]] == 0] <- 1.0e-323
      Q_old[[j]][Q_old[[j]] == 0] <- 1.0e-323
    }
    pi_2_old <- rdirichlet(alpha_2_old/M_max + table_int(unlist(z_sl_old), M_max))
    
    
    # (4) Update the hyperparameter for the transition probabilities
    for (j in 1:p){
      alpha_temp <- exp(rnorm(1, log(alpha_old[j]), sd_MH_alpha[j]))
      lu <- log(runif(1))
      acc_rate <- min(c(0, (d_j[j] + 1) * (lgamma(alpha_temp) - d_j[j] * lgamma(alpha_temp/d_j[j])) + 
                          (alpha_temp/d_j[j] - 1) * (sum(log(Q_old[[j]])) + sum(log(pi_0_old[[j]]))) +
                          dgamma(alpha_temp, hypers$a_alpha, hypers$b_alpha, log = T) +
                          log(alpha_temp) -
                          (d_j[j] + 1) * (lgamma(alpha_old) - d_j[j] * lgamma(alpha_old[j]/d_j[j])) - 
                          (alpha_old[j]/d_j[j] - 1) * (sum(log(Q_old[[j]])) + sum(log(pi_0_old[[j]]))) -
                          dgamma(alpha_old[j], hypers$a_alpha, hypers$b_alpha, log = T) -
                          log(alpha_old[j])))
      if (lu < acc_rate){
        alpha_old[j] <- alpha_temp
        acc_alpha[j] <- acc_alpha[j] + 1
      }
    }
    
    
    # (5) Update the cluster specific smoothness parameters
    sigma2_1_temp <- exp(rnorm(1, log(sigma2_1_old), sd_MH_sigma))
    lu <- log(runif(1))
    acc_rate <- min(c(0, -0.5 * n_RSS * log(sigma2_1_temp) -
                        0.5/sigma2_1_temp * RSS +
                        dhalfcauhy(sqrt(sigma2_1_temp), hypers$sigma_s, log = T) +
                        log(sigma2_1_temp) +
                        0.5 * n_RSS * log(sigma2_1_old) +
                        0.5/sigma2_1_old * RSS -
                        dhalfcauhy(sqrt(sigma2_1_old), hypers$sigma_s, log = T) -
                        log(sigma2_1_old)))
    if (lu < acc_rate){
      sigma2_1_old <- sigma2_1_temp
      acc_sigma <- acc_sigma + 1
    }
    
    
    # Compute residuals: data minus main effects	
    res_u <- y - B_beta
    
    # (6) Update random effects	
    reff <- sample_reff(res_u, B, P_u, ind, sigma2_y_old, sigma2_us_old, 	
                        sigma2_ua_old)	
    beta_ui_old <- reff$beta_ui	
    B_beta_u <- as.numeric(reff$B_beta_u)	
    
    
    # (7) Update random effects variances	
    sigma2_ua_prop <- exp(rnorm(1, mean = log(sigma2_ua_old), sd = sd_MH_sigma_ua))	
    acc_rate <- min(c(0, 0.5 * n_ind * log(det(diag(K)/sigma2_ua_prop + P_u/sigma2_us_old)) -	
                        0.5 * sum(diag(beta_ui_old %*% tcrossprod(diag(K)/sigma2_ua_prop, beta_ui_old))) +	
                        dhalfcauhy(sqrt(sigma2_ua_prop), hypers$sigma_s, TRUE) + 	
                        log(sigma2_ua_prop) -	
                        0.5 * n_ind * log(det(diag(K)/sigma2_ua_old + P_u/sigma2_us_old)) +	
                        0.5 * sum(diag(beta_ui_old %*% tcrossprod(diag(K)/sigma2_ua_old, beta_ui_old))) -	
                        dhalfcauhy(sqrt(sigma2_ua_old), hypers$sigma_s, TRUE) - 	
                        log(sigma2_ua_old)))	
    lu <- log(runif(1))	
    if (lu < acc_rate){	
      sigma2_ua_old <- sigma2_ua_prop	
      acc_sigma_ua <- acc_sigma_ua + 1	
    }	
    sigma2_us_prop <- exp(rnorm(1, mean = log(sigma2_us_old), sd = sd_MH_sigma_us))	
    acc_rate <- min(c(0, 0.5 * n_ind * log(det(diag(K)/sigma2_ua_old + P_u/sigma2_us_prop)) -	
                        0.5 * sum(diag(beta_ui_old %*% tcrossprod(P_u/sigma2_us_prop, beta_ui_old))) +	
                        dhalfcauhy(sqrt(sigma2_us_prop), hypers$sigma_s, TRUE) + 	
                        log(sigma2_us_prop) -	
                        0.5 * n_ind * log(det(diag(K)/sigma2_ua_old + P_u/sigma2_us_old)) +	
                        0.5 * sum(diag(beta_ui_old %*% tcrossprod(P_u/sigma2_us_old, beta_ui_old))) -	
                        dhalfcauhy(sqrt(sigma2_us_old), hypers$sigma_s, TRUE) - 	
                        log(sigma2_us_old)))	
    if (lu < acc_rate){	
      sigma2_us_old <- sigma2_us_prop	
      acc_sigma_us <- acc_sigma_us + 1	
    }
    
    
    # (8) Update the variance of the data (global)
    RSS <- as.numeric(crossprod(y - B_beta - B_beta_u))
    sigma2_y_old <- 1/rgamma(1, hypers$a_sigma + n/2, hypers$b_sigma + 0.5 * RSS)
    
    
    # (9) Update the dimensionality penalty phi_old
    for (j in 1:p){
      phi_prop <- exp(rnorm(1, log(phi_old[j]), sd_MH_phi[j]))
      lu <- log(runif(1))
      acc_rate <- min(c(0, - phi_prop * sum(l_k_old[j,]) - 
                          K * log(sum(exp(- phi_prop * 1:d_j[j]))) + 
                          dgamma(phi_prop, hypers$a_phi, hypers$b_phi, log = T) +
                          log(phi_prop) - (
                            - phi_old[j] * sum(l_k_old[j,]) - 
                              K * log(sum(exp(- phi_old[j] * 1:d_j[j]))) + 
                              dgamma(phi_old[j], hypers$a_phi, hypers$b_phi, log = T) +
                              log(phi_old[j]))))
      if (lu < acc_rate){
        phi_old[j] <- phi_prop
        acc_phi[j] <- acc_phi[j] + 1
      }
    }
    
    
    # Save parameters in the chain
    if ( (iter > burnin) && (iter %% thin == 0) ){
      for (i in 1:n_ind){
        post_reff_gr[,i,it] <- Bgrid %*% beta_ui_old[i,]
      }
      phi_chain[it,] <- phi_old
      sigma2_1_chain[it] <- sigma2_1_old
      l_k_chain[it,,] <- l_k_old
      post_data[it,] <- B_beta + B_beta_u
      post_pred[it,,] <- B_beta_pred
      post_mean[it,] <- B_beta
      sigma2_1_chain[it] <- sigma2_1_old
      sigma2_y_chain[it] <- sigma2_y_old
      alpha_chain[it,] <- alpha_old
      RSS_chain[it] <- RSS
      sigma2_us_chain[it] <- sigma2_us_old	
      sigma2_ua_chain[it] <- sigma2_ua_old
      
      it <- it + 1
    }
    setTxtProgressBar(pb, iter/Niter)
  }
  
  return(list('post_data' = post_data, 'post_mean' = post_mean, 
              'l_k' = l_k_chain, 'sigma2_1' = sigma2_1_chain, 
              'sigma2_y' = sigma2_y_chain, 'alpha' = alpha_chain, 
              'sigma2_us_chain' = sigma2_us_chain, 
              'sigma2_ua_chain' = sigma2_ua_chain, 
              'post_pred' = post_pred, 
              'RSS' = RSS_chain, 'post_reff_gr' = post_reff_gr))
  
}


