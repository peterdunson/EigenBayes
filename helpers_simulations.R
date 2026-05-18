

source('competitors/FACTOR_ANALYSIS/FACTOR_CODE_update.R')
source('fable.R')
library(FABLE)
Rcpp::sourceCpp("competitors/FABLE/src/updated-FABLE-functions.cpp")
library(infinitefactor)



library(LaplacesDemon)
library(MASS)
library(matrixStats)
library(readxl)
library(truncnorm)

source('eigenbayes_functions.R')



gen_factor_data <- function(n = 300, p = 200, r=10,
                            psi_min = 0.5, psi_max = 1.5,
                            loading_scales = rep(1, 10), adapt_to_outcome = F, sparse=F) {
  # structured Lambda (p x r)
  Lambda <- matrix(0, p, r)
  idx <- (1:p) / p
  for (j in 1:r) {
    base <- exp(-3 * abs(idx - j/(r + 1L)))
    Lambda[, j] <- loading_scales[1] * (base + 0.25 * rnorm(p))
  }
  # normalize columns to have ||col||^2 ≈ p
  #Lambda <- scale(Lambda, center = FALSE, scale = sqrt(colSums(Lambda^2) / p))
  if(! adapt_to_outcome){
    Lambda <- matrix(rnorm(p*r), ncol=r) %*% diag(as.vector(loading_scales))
  }
  #Lambda <- cbind(Lambda, loading_scale_2 * matrix(rnorm(p*r_2), ncol=r_2))
  
  if(adapt_to_outcome){
    Lambda <- matrix(NA, p, r)
    i_1 <- floor(p / 3)
    i_2 <- floor(p * 2 / 3)
    Lambda[1:i_1,] <- sqrt(3) * matrix(rnorm(i_1*r), ncol=r) %*% diag(as.vector(loading_scales))
    Lambda[(i_1+1):i_2,] <- matrix(rnorm((i_2 - i_1)*r), ncol=r) %*% diag(as.vector(loading_scales))
    Lambda[-c(1:i_2),] <- 1/sqrt(3) * matrix(rnorm((p - i_2)*r), ncol=r) %*% diag(as.vector(loading_scales))
  }
  
  if(sparse){
    sparse_idx <- sample(1:p, floor(p/3))
    Lambda[sparse_idx,] <- 0
  }else{
    sparse_idx <- c()
  }
  
  # factors F ~ N(0, I_r) across rows
  Fmat <- matrix(rnorm(n * (r)), n, r)
  
  # heteroskedastic idiosyncratic variances per column (variable)
  psi <- runif(p, psi_min, psi_max)
  E <- matrix(rnorm(n * p), n, p)
  E <- sweep(E, 2, sqrt(psi), "*")  # column-wise variance
  
  X <- Fmat %*% t(Lambda) + E
  list(Y = X, Lambda = Lambda, Sigma = psi, F_mat = Fmat, sparse_idx=sparse_idx)
}



gen_factor_data_2 <- function(n = 300, p = 200, r_1 = 3,r_2 = 3,
                              psi_min = 0.5, psi_max = 1.5,
                              loading_scale_1 = 1.0, loading_scale_2 = 0.1,
                              adapt_to_outcome = F) {
  # structured Lambda (p x r)
  Lambda <- matrix(0, p, r_1 +r_2)
  idx <- (1:p) / p
  Lambda <- loading_scale_1 * matrix(rnorm(p*r_1), ncol=r_1)
  Lambda <- cbind(Lambda, loading_scale_2 * matrix(rnorm(p*r_2), ncol=r_2))
  if(adapt_to_outcome){
    i_1 <- floor(p / 3)
    i_2 <- floor(p * 2 / 3)
    Lambda[1:i_1,] <- sqrt(3) *  Lambda[1:i_1,]
    Lambda[-c(1:i_2),] <- 1/sqrt(3) * Lambda[-c(1:i_2),]
  }
  # factors F ~ N(0, I_r) across rows
  Fmat <- matrix(rnorm(n * (r_1 + r_2)), n, r_1 + r_2)
  
  # heteroskedastic idiosyncratic variances per column (variable)
  psi <- runif(p, psi_min, psi_max)
  E <- matrix(rnorm(n * p), n, p)
  E <- sweep(E, 2, sqrt(psi), "*")  # column-wise variance
  
  X <- Fmat %*% t(Lambda) + E
  list(Y = X, Lambda = Lambda, Sigma = psi, F_mat = Fmat)
}


barigozzi_cho_est <- function(Y, k, c_w=1.1){
  n <- nrow(Y); p <- ncol(Y)
  s_Y <- svd(Y, nu=k, nv=k)
  V <- s_Y$v[, 1:k, drop=F]
  V_scaled <- V
  c_w <- c_w * max(abs(V[,1]))
  for(l in 1:k){
    v_l <- max(1, 1 / c_w * max(abs(V[,l])))
    V_scaled[,l] <-  V_scaled[,l] / v_l
  }
  Y_hat <- Y %*%  tcrossprod(V_scaled)
  s_Y <- svd(Y_hat, nu=k, nv=k)
  M_hat <- s_Y$u[,1:k, drop=F] * sqrt(n)
  Lambda_hat <- s_Y$v[,1:k, drop=F] %*% diag(as.vector(s_Y$d[1:k, drop=F])) / sqrt(n)
  sigmas_hat_2 <- colSums( (Y - tcrossprod(M_hat, Lambda_hat))^2) / (n - k)
  sigmas_hat <- colSums( (Y - tcrossprod(M_hat, Lambda_hat))^2) / n
  
  return(list(
    Y_hat = Y_hat,
    M_hat = M_hat,
    Lambda_hat = Lambda_hat,
    sigmas_sq_hat = sigmas_hat,
    sigmas_sq_hat_2 = sigmas_hat_2
    
  ))
}


eigenvals_cov <- function(Y, center = FALSE) {
  Y <- as.matrix(Y)
  if (center) Y <- scale(Y, center = TRUE, scale = FALSE)
  n <- nrow(Y)
  #S <- crossprod(Y) / n
  #as.numeric(eigen(S, symmetric = TRUE, only.values = TRUE)$values)
  return(svd(Y)$d^2 / n)
}

estimate_r_bn <- function(Y, rmax = NULL, g = NULL, center = FALSE) {
  Y <- as.matrix(Y)
  n <- nrow(Y)
  p <- ncol(Y)
  
  mu <- eigenvals_cov(Y, center = center)
  
  if (is.null(rmax)) rmax <- floor(sqrt(min(n, p)))
  rmax <- max(1L, min(as.integer(rmax), p - 1L))
  
  if (is.null(g)) {
    g <- (p + n) * log(min(p, n)) / (p * n)
  }
  
  IC <- numeric(rmax)
  for (q in 1:rmax) {
    l_max <- min(n, p)
    
    tail_mean <- sum(mu[(q + 1):l_max]) / p
    tail_mean <- max(tail_mean, .Machine$double.eps)
    IC[q] <- log(tail_mean) + q * g
  }
  
  list(
    k_hat = which.min(IC),
    IC = IC,
    eigenvalues = mu,
    rmax = rmax,
    g = g
  )
}

estimate_r_ah <- function(Y, rmax = NULL, center = FALSE) {
  Y <- as.matrix(Y)
  n <- nrow(Y)
  p <- ncol(Y)
  
  mu <- eigenvals_cov(Y, center = center)
  
  if (is.null(rmax)) rmax <- floor(sqrt(min(p, n)))
  rmax <- max(1L, min(as.integer(rmax), p - 2L))
  l_max <- min(n, p)
  tail_sum <- sapply(1:(l_max - 1), function(q) sum(mu[(q + 1):l_max]))
  tail_sum <- pmax(tail_sum, .Machine$double.eps)
  
  mu_star <- mu[1:(l_max - 1)] / tail_sum
  
  GR <- numeric(rmax)
  for (q in 1:rmax) {
    num <- log(1 + mu_star[q])
    den <- log(1 + mu_star[q + 1])
    den <- max(den, .Machine$double.eps)
    GR[q] <- num / den
  }
  
  list(
    k_hat = which.max(GR),
    GR = GR,
    mu_star = mu_star,
    eigenvalues = mu,
    rmax = rmax
  )
}

compute_metrics <- function(est, compute_coverage_=F, idx_cvg=1:100, idx_fr=1:100){
  res <- c()
  est_cov <- tcrossprod(est$Lambda_hat)
  
  res[1] <- fro_rel_err(est_cov, Lambda_outer_0)
  res[2] <- fro_rel_err(est_cov[idx_fr, idx_fr], Lambda_outer_0[idx_fr, idx_fr])
  
  est_cov  <- est_cov + diag(est$sigmas_sq_hat)
  res[3] <- fro_rel_err(est_cov, Theta_0)
  res[4] <- fro_rel_err(est_cov[idx_fr, idx_fr], Theta_0[idx_fr, idx_fr])
  res[5] <- fro_rel_err(est$Y_hat, X)
  res[6] <- fro_rel_err(est$Y_hat[,idx_fr], X[,idx_fr])
  
  if(compute_coverage_){
    eb_cis <- eigenbayes_approx_ci(est, subsample_index=idx_cvg)
    cov <- compute_coverage(Theta_0[idx_cvg,idx_cvg], eb_cis, subsample_index=1:length(idx_cvg))
    res[7] <- mean(cov$coverage_confidence_intervals)
    res[8] <- mean(cov$length_confidence_intervals)
    res[9] <- mean(cov$coverage_credible_intervals)
    res[10] <- mean(cov$length_credible_intervals)
  }
  else{
    res = c(res, rep(0,4))
  }
  
  return(res)
}

compute_metrics_fable <- function(fit, compute_coverage_=F, idx_cvg=1:100, idx_fr=1:100){
  res <- c()
  #res[1] <- fro_rel_err(fit$Lambda_outer, Lambda_outer_0)
  #res[2] <- fro_rel_err(fit$Lambda_outer[idx_fr, idx_fr], Lambda_outer_0[idx_fr, idx_fr])
  
  res[3] <- fro_rel_err(fit$FABLEPostMean, Theta_0)
  res[4] <- fro_rel_err(fit$FABLEPostMean[idx_fr, idx_fr], Theta_0[idx_fr, idx_fr])
  
  #res[5] <- fro_rel_err(fit$Y_hat, X)
  #res[6] <- fro_rel_err(fit$Y_hat[, idx_fr], X[, idx_fr])
  
  if(compute_coverage_){
    FABLESamples = FABLEPosteriorSampler(data$Y, gamma0 = 1, delta0sq = 1, maxProp = 0.95, MC = 1000)
    #fable_cov_samples <- construct_fable_cov_samples(FABLESamples, idx_cvg)
    #fable_cis <- fable_monte_carlo_ci(fable_cov_samples)
    ptm <- proc.time()
    svdmod = svd(data$Y)
    U_Y = svdmod$u
    V_Y = svdmod$v
    svalsY = svdmod$d
    kEst = CPPRankEstimator(data$Y, U_Y, V_Y, svalsY, 50)
    kEst
    FABLEHypPars = FABLEHyperParameters(data$Y, U_Y, V_Y, svalsY, kEst)
    covCorrectEntries = CPPcov_correct_matrix(FABLEHypPars$SigmaSqEstimate,
                                              FABLEHypPars$G)
    varInflation = mean(covCorrectEntries)^2
    varInflation
    CPPSamplingOutput = CPPFABLESampler(data$Y, 
                                        1, 
                                        1, 
                                        1000,
                                        U_Y,
                                        V_Y,
                                        svalsY,
                                        kEst,
                                        varInflation)
    
    CPPPostProcess = CCFABLEPostProcessingSubmatrix(CPPSamplingOutput,
                                                    0.05,
                                                    idx_cvg)
    
    
    time_uq <- proc.time() - ptm; 
    fable_cis <- list()
    p <- ncol(data$Y)
    fable_cis$credible_intervals <- array(NA, dim=c(2, p, p) )
    fable_cis$credible_intervals[1,,] <- CPPPostProcess$LowerQuantileMatrix
    fable_cis$credible_intervals[2,,] <- CPPPostProcess$UpperQuantileMatrix
    
    
    fable_cov <- compute_coverage(Theta_0[idx_cvg,idx_cvg], fable_cis, confidence_intervals=F, subsample_index=1:length(idx_cvg))
    
    res[1] <- fro_rel_err(fit$Lambda_outer, Lambda_outer_0)
    res[2] <- fro_rel_err(fit$Lambda_outer[idx_fr, idx_fr], Lambda_outer_0[idx_fr, idx_fr])
    
    Y_hat <- fable_low_rank_signal(data$Y, FABLESamples)
    
    
    res[5] <- fro_rel_err(Y_hat, X)
    res[6] <- fro_rel_err(Y_hat[, idx_fr], X[, idx_fr])
    
    res[7] <- 0
    res[8] <- 0
    res[9] <- mean(fable_cov$coverage_credible_intervals)
    res[10] <- mean(fable_cov$length_credible_intervals)
  }
  else{
    res = c(res, rep(0,4))
  }
  return(res)
}

.fro <- function(M) sqrt(sum(M * M))
fro_rel_err <- function(Ahat, Atrue) {
  num <- .fro(Ahat - Atrue); den <- .fro(Atrue)
  if (den == 0) return(NA_real_)
  num / den
}


rotate_est <- function(Y, k, lambda1=0.001){
  p <- ncol(Y)
  n <- nrow(Y)
  startB <- matrix(rnorm(p*k), p, k)
  alpha <- 1/p
  lambda1 <- 0.001
  epsilon <- 0.05
  lambda0=5
  start <- list(B=startB, sigma=rep(1,k), theta=rep(0.5,k))
  rotate_fit <- FACTOR_ROTATE(
    Y, lambda0, lambda1, start, k, epsilon, alpha,TRUE,TRUE,100,T, plot=F
  )
  #rotate_fit$Lambda_outer <- tcrossprod(rotate_fit$B)
  rotate_fit$Lambda_hat <- rotate_fit$B
  rotate_fit$sigmas_sq_hat <- rotate_fit$sigma^2
  
  rotate_fit$M_hat <- Y %*% diag(1/rotate_fit$sigmas_sq_hat) %*%  rotate_fit$Lambda_hat %*%
    solve(t(rotate_fit$Lambda_hat) %*% diag(1/rotate_fit$sigmas_sq_hat) %*%  rotate_fit$Lambda_hat + diag(1,k,k))
  rotate_fit$Y_hat <- rotate_fit$M_hat %*% t(rotate_fit$Lambda_hat)
  
  return(rotate_fit)
}

spectral_est <- function(Y, s_Y, k){
  U <- s_Y$u[,1:k, drop=F]
  V <- s_Y$v[,1:k, drop=F]
  D <- diag(as.vector(s_Y$d[1:k, drop=F]))
  sigmas_sq_2 <- colSums((Y - tcrossprod(U) %*% Y)^2 ) /  (n - k)
  sigmas_sq <- colSums((Y - tcrossprod(U) %*% Y)^2 ) /  n
  Lambda_hat <-  V %*% D / sqrt(n)
  return(list(Lambda_hat = Lambda_hat, 
              sigmas_sq_hat = sigmas_sq,
              sigmas_sq_hat_2 = sigmas_sq_2,
              Y_hat = U %*% D %*% t(V)))
}

