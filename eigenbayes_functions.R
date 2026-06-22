library(Rcpp)
library(RcppArmadillo)
sourceCpp('helper_functions_eigenbayes.cpp')


fit_eigenbayes <- function(Y, H,v_0=5){
  s_Y <- svd(Y)
  res <- eigenbayes_point_est(Y, s_Y, H, v_0)
  return(res)
}

eigenbayes_point_est <- function(
    Y, s_Y, k,v_0=5,
){
  p <- ncol(Y)
  n <- nrow(Y)
  U <- s_Y$u[,1:k]
  V <- s_Y$v[,1:k]
  d <- s_Y$d[1:k]
  n <- nrow(Y)
  M <- U * sqrt(n)
  d_bar_sq <- mean(s_Y$d[-c(1:k)]^2)
  d_bar_sq <- s_Y$d[k+1]^2
  D <- diag(as.vector(d))
  sigmas_sq <- colSums((Y - tcrossprod(U) %*% Y)^2 ) /  (n - k)
  sigma_sq_0 <- mean(sigmas_sq)
  signals_magnitude <- (d^2 - d_bar_sq)^alpha
  phis_sq <-(signals_magnitude) / sum(signals_magnitude)  * k
  
  taus_sq <- diag(V %*% diag(d^2 - d_bar_sq) %*% t(V)) /
    (n * k * sigmas_sq) 
  Lambda_hat <- compute_point_estimate_Lambda_hat(V, D, n, taus_sq, phis_sq)
  Y_hat <- tcrossprod(M, Lambda_hat)
  sigmas_sq_mean <- (v_0 *sigma_sq_0 + colSums((Y - Y_hat)^2)) / (n + v_0 -2)
  
  return(list(Lambda_hat = Lambda_hat,
              sigmas_sq_mean = sigmas_sq_mean,
              sigmas_sq_hat = sigmas_sq,
              sigma_sq_0 =sigma_sq_0,
              taus_sq = taus_sq,
              phis_sq = phis_sq,
              M_hat = M,
              Y_hat = Y_hat,
              v_n = v_0 + n)
  )
}


eigenbayes_covariance_est <- function(
    eb_fit, only_low_rank=F, posterior_mean=T, spherical_idyosincratic_covariance=T
){
  Lambda_outer <- tcrossprod(eb_fit$Lambda_hat)
  p <- ncol(Lambda_outer)
  if(only_low_rank){
    if(! posterior_mean){
      return(Lambda_outer)
    }
  }
  if(is.null(eb_fit$rho_sq)){
    eb_fit <- compute_rho_sq(eb_fit, 1:p, spherical_idyosincratic_covariance=spherical_idyosincratic_covariance)
  }
  correction <- c()
  for(j in 1:p){
    correction[j] <- eb_fit$sigmas_sq_mean * sum(1/(n + 1/eb_fit$taus_sq[j]*1/eb_fit$phis_sq)) * eb_fit$rho_sq
  }
  Lambda_outer <- Lambda_outer + diag(correction) 
  if(only_low_rank){
    return(Lambda_outer)
  }
  covariance <- Lambda_outer + diag(eb_fit$sigmas_sq_mean)
  return(list(Lambda_outer=Lambda_outer, covariance=covariance))
}


eigenbayes_posterior_samples <- function(
    eb_fit, n_MC=500, subsample_index=NULL, spherical_idyosincratic_covariance=T
){
  p <- nrow(eb_fit$Lambda_hat); k <- ncol(eb_fit$Lambda_hat)
  if(is.null(subsample_index)){
    subsample_index <- 1:p
  }
  p_index <- length(subsample_index)
  sigmas_sq_samples <- sapply(
    eb_fit$sigmas_sq_mean[subsample_index], function(x) (rgamma(n_MC, eb_fit$v_n, x))
  ) # n_MC * p
  sigmas_sq_samples <- 1/sigmas_sq_samples
  if(is.null(eb_fit$rho_sq)){
    eb_fit <- compute_rho_sq(eb_fit, subsample_index, spherical_idyosincratic_covariance=spherical_idyosincratic_covariance)
  }
  #Lambda_samples <- array(NA, dim=c(n_MC, p_index, k))
  #for(t in 1:n_MC){
  #  for(j in 1:p_index){
  #    lambda_j <- eb_fit$Lambda_hat[subsample_index[j],] + rnorm(k, 0, 1) * sqrt(eb_fit$rho_sq * sigmas_sq_samples[t,j]) * 1/sqrt(n + 1/eb_fit$taus_sq[subsample_index[j]]*1/eb_fit$phis_sq)
  #    Lambda_samples[t, j, ] <- lambda_j
  #   }
  #}
  Lambda_samples <- sample_Lambda(
    eb_fit$Lambda_hat, subsample_index, eb_fit$taus_sq,eb_fit$eb_fit$phis_sq, 
    sigmas_sq_samples, eb_fit$rho_sq, n
  )
  
  return(list(sigmas_sq_samples=sigmas_sq_samples, Lambda_samples=Lambda_samples))
}



compute_rho_sq <- function(eb_fit, subsample_index=NULL, spherical_idyosincratic_covariance=T){
  
  p <- nrow(eb_fit$Lambda_hat); k <- ncol(eb_fit$Lambda_hat)
  if(is.null(subsample_index)){
    subsample_index <- 1:p
  }
  p_index <- length(subsample_index)
  
  if(is.null(eb_fit$Lambda_outer)){
    eb_fit$Lambda_outer <-  tcrossprod(eb_fit$Lambda_hat[,])
  }
  
  if(spherical_idyosincratic_covariance){
    sigmas_sq = rep(mean(eb_fit$sigmas_sq_mean), p)
  } else {
    sigmas_sq = eb_fit$sigmas_sq_mean
  }
  
  correction <- compute_correction(
    eb_fit$Lambda_outer, sigmas_sq, subsample_index)
  
  rho <- mean(correction[upper.tri(correction, diag=T)])
  eb_fit$rho_sq <- rho^2
  
  return(eb_fit)
}



eigenbayes_approx_ci <- function(
    eb_fit, alpha=0.05, return_confidence_intervals=T, return_credible_intervals=T,
    spherical_idyosincratic_covariance=T, subsample_index=NA){
  
  if(any(is.na(subsample_index))){
    p <- nrow(eb_fit$Lambda_hat)
    subsample_index <- 1:p
  }
  
  
  n <- nrow(eb_fit$M_hat)
  p <- nrow(eb_fit$Lambda_hat)
  if(is.null(eb_fit$rho_sq)){
    eb_fit <- compute_rho_sq(eb_fit, subsample_index, spherical_idyosincratic_covariance=spherical_idyosincratic_covariance)
  }
  
  output <- list()
  point_estimator <- eb_fit$Lambda_outer[subsample_index, subsample_index] + diag(eb_fit$sigmas_sq_mean[subsample_index])
  output$covariance_estimate  <- point_estimator
  
  if(spherical_idyosincratic_covariance){
    sigmas_sq = rep(mean(eb_fit$sigmas_sq_mean), p)
  } else {
    sigmas_sq = eb_fit$sigmas_sq_mean
  }
  
  p <- length(subsample_index)
  
  if(return_confidence_intervals){
    confidence_intervals <- array(NA, dim=c(2, p, p))
    sds_clt <- compute_sds_clt(eb_fit$Lambda_outer[subsample_index, subsample_index], sigmas_sq[subsample_index])
    dev <- sds_clt * qnorm(1-alpha/2) / sqrt(n)
    confidence_intervals[1,,] = point_estimator - dev
    confidence_intervals[2,,] = point_estimator + dev
    output$confidence_intervals <- confidence_intervals 
  }
  
  if(return_credible_intervals){
    credible_intervals <- array(NA, dim=c(2, p, p))
    sds_bvm <- compute_sds_bvm(eb_fit$Lambda_outer[subsample_index,subsample_index], sigmas_sq[subsample_index], eb_fit$rho_sq)
    dev <- sds_bvm * qnorm(1-alpha/2)/ sqrt(n)
    credible_intervals[1,,] = point_estimator - dev
    credible_intervals[2,,] = point_estimator + dev
    output$credible_intervals <- credible_intervals 
  }
  
  return(output)
}

compute_coverage <- function(
    Theta_0, cis, confidence_intervals=T, credible_intervals=T, subsample_index=NA){
  
  if(any(is.na(subsample_index))){
    p <- ncol(Theta_0)
    subsample_index <- 1:p
  }
  
  idx <- upper.tri(Theta_0[subsample_index, subsample_index], diag=T)
  res <- list()
  
  if(confidence_intervals){
    res$coverage_confidence_intervals <- (Theta_0[idx] > (cis$confidence_intervals[1,subsample_index, subsample_index][idx])) & (Theta_0[idx] < (cis$confidence_intervals[2,subsample_index, subsample_index][idx]))
    res$length_confidence_intervals <- cis$confidence_intervals[2,subsample_index, subsample_index][idx] - cis$confidence_intervals[1,subsample_index, subsample_index][idx]
    print(paste0('mean coverage conf int : ', mean(res$coverage_confidence_intervals)))
    print(paste0('mean length conf int : ', mean(res$length_confidence_intervals)))
    
  }
  
  if(credible_intervals){
    res$coverage_credible_intervals <- (Theta_0[idx] >  (cis$credible_intervals[1,subsample_index, subsample_index][idx])) & (Theta_0[idx] < (cis$credible_intervals[2,,][idx]))
    res$length_credible_intervals <- cis$credible_intervals[2,subsample_index, subsample_index][idx] - cis$credible_intervals[1,subsample_index, subsample_index][idx]
    print(paste0('mean coverage cred int : ', mean(res$coverage_credible_intervals)))
    print(paste0('mean length cred int : ', mean(res$length_credible_intervals)))
    
  }
  
  return(res)
}


latent_factor_mean_old <- function(Y, Lambda_samples, sigmas_sq_samples){
  k <- dim(Lambda_samples)[3]
  p <- dim(Lambda_samples)[2]
  n_MC <- dim(Lambda_samples)[1]
  n <- nrow(Y)
  
  latent_factors_conditional_mean <- array(NA, dim=c(n_MC, n, k))
  for(t in 1:n_MC){
    variance_t <- t(Lambda_samples[t,,]) %*% diag(1/sigmas_sq_samples[t,]) %*% Lambda_samples[t,,] + diag(1, k, k)
    variance_t <- solve(variance_t)
    latent_factors_conditional_mean[t,,] <- latent_factors_full_conditional_mean(
      Y, Lambda_samples[t,,], sigmas_sq[t,], variance_t
    )
  }
  return(apply(latent_factors_conditional_mean, c(2,3), mean))
  
}


latent_factor_full_conditional_old <- function(mean_lf, variance_lf){
  n <- nrow(mean_lf)
  k <- ncol(mean_lf)
  return(rmvnorm(n, rep(0, k), variance_lf) + mean_lf)
}




