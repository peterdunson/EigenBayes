

source('competitors/FACTOR_ANALYSIS/FACTOR_CODE_update.R')
source('fable.R')
library(infinitefactor)



library(LaplacesDemon)
library(MASS)
library(matrixStats)
library(readxl)
library(truncnorm)

library(Rcpp)
library(RcppArmadillo)
sourceCpp('helper_functions_eigenbayes.cpp')


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
  c_w <- c_w * max(V[,1])
  for(l in 1:k){
    v_l <- max(1, 1 / c_w * max(V[,l]))
    V_scaled[,l] <-  V_scaled[,l] / v_l
  }
  Y_hat <- Y %*%  tcrossprod(V_scaled)
  s_Y <- svd(Y_hat, nu=k, nv=k)
  M_hat <- s_Y$u[,1:k, drop=F] * sqrt(n)
  Lambda_hat <- s_Y$v[,1:k, drop=F] %*% diag(as.vector(s_Y$d[1:k, drop=F])) / sqrt(n)
  sigmas_hat <- colSums( (Y - tcrossprod(M_hat, Lambda_hat))^2) / (n - k)
  return(list(
    Y_hat = Y_hat,
    M_hat = M_hat,
    Lambda_hat = Lambda_hat,
    sigmas_sq_hat = sigmas_hat
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

compute_metrics <- function(est, compute_coverage=F, idx_cvg=1:100, idx_fr=1:100){
  res <- c()
  est_cov <- tcrossprod(est$Lambda_hat)
  
  res[1] <- fro_rel_err(est_cov, Lambda_outer_0)
  res[2] <- fro_rel_err(est_cov[idx_fr, idx_fr], Lambda_outer_0[idx_fr, idx_fr])
  
  est_cov  <- est_cov + diag(est$sigmas_sq_hat)
  res[3] <- fro_rel_err(est_cov, Theta_0)
  res[4] <- fro_rel_err(est_cov[idx_fr, idx_fr], Theta_0[idx_fr, idx_fr])
  res[5] <- fro_rel_err(est$Y_hat, X)
  res[6] <- fro_rel_err(est$Y_hat[,idx_fr], X[,idx_fr])
  
  if(compute_coverage){
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

compute_metrics_fable <- function(fit, compute_coverage=F, idx_cvg=1:100, idx_fr=1:100){
  res <- c()
  #res[1] <- fro_rel_err(fit$Lambda_outer, Lambda_outer_0)
  #res[2] <- fro_rel_err(fit$Lambda_outer[idx_fr, idx_fr], Lambda_outer_0[idx_fr, idx_fr])
  
  res[3] <- fro_rel_err(fit$FABLEPostMean, Theta_0)
  res[4] <- fro_rel_err(fit$FABLEPostMean[idx_fr, idx_fr], Theta_0[idx_fr, idx_fr])

  #res[5] <- fro_rel_err(fit$Y_hat, X)
  #res[6] <- fro_rel_err(fit$Y_hat[, idx_fr], X[, idx_fr])

  
  
  if(compute_coverage){
    FABLESamples = FABLEPosteriorSampler(data$Y, gamma0 = 1, delta0sq = 1, maxProp = 0.95, MC = 1000)
    fable_cov_samples <- construct_fable_cov_samples(FABLESamples, idx_cvg)
    fable_cis <- fable_monte_carlo_ci(fable_cov_samples)
    fable_cov <- compute_coverage(Theta_0[idx_cvg,idx_cvg], fable_cis, confidence_intervals=F, subsample_index=1:length(idx_cvg))
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
    Y, lambda0, lambda1, start, k, epsilon, alpha,TRUE,TRUE,100,F
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
  sigmas_sq <- colSums((Y - tcrossprod(U) %*% Y)^2 ) /  (n - k)
  Lambda_hat <-  V %*% D / sqrt(n)
  return(list(Lambda_hat = Lambda_hat, 
             sigma_hat = sigmas_sq,
              X = U %*% D %*% t(V)))
}



n_sim <- 50
scenario <- 1
if(scenario == 1 | scenario == 2 ){
  k <- 10
  p <- 1000
  n <- 500
  if(scenario == 2){
    n <- 1000
  }
}
if(scenario == 3 | scenario == 4 ){
  k <- 10
  p <- 5000
  n <- 500
  if(scenario == 4){
    n <- 1000
  }
}

subsample_index<- 1:100 


bc_true_k <- data.frame()
bc_k_bn <- data.frame()
bc_k_ah <- data.frame()
bc_k_jic <- data.frame()
bc_k_jic_over <- data.frame()
bc_k_over <- data.frame()

eb_true_k <- data.frame()
eb_k_bn <- data.frame()
eb_k_ah <- data.frame()
eb_k_jic <- data.frame()
eb_k_jic_over <- data.frame()
eb_k_over <- data.frame()

rotate_true_k <- data.frame()
rotate_k_bn <- data.frame()
rotate_k_ah <- data.frame()
rotate_k_jic <- data.frame()
rotate_k_jic_over <- data.frame()
rotate_k_over <- data.frame()

pca_true_k <- data.frame()
pca_k_bn <- data.frame()
pca_k_ah <- data.frame()
pca_k_jic <- data.frame()
pca_k_jic_over <- data.frame()
pca_k_over <- data.frame()

fable_k_hat <- data.frame()

test_barigozzi_cho <- T; test_eigenbayes <- T; test_fable <- T; test_pca <- T; test_rotate <- T

adapt_to_outcome <- T; sparse <- F

for(sim in 1:n_sim){
  set.seed(sim)
  print(sim)
  data <- gen_factor_data(
    n=n, p=p, psi_min=1, psi_max=5, r=k, 
    loading_scales=seq(0.15, 0.25, length.out = k),
    adapt_to_outcome = adapt_to_outcome, sparse = sparse
  )
  
  
  idx_cvg=sample(1:p, 100)
  idx_fr=(floor(p * 2 / 3)+1):p
  if(sparse){
    idx_fr = data$sparse_idx 
  }
  
  #compute_coverage=T, idx_cvg=idx_cvg, idx_fr=idx_fr
  
  
  Lambda_outer_0 <- tcrossprod(data$Lambda)
  Theta_0 <- Lambda_outer_0 + diag(data$Sigma)
  X <- data$F_mat %*% t(data$Lambda)
  
  k_bn <- estimate_r_bn(data$Y)
  k_ah <- estimate_r_ah(data$Y)
  k_over <- sqrt(n)
  
  if(test_fable){
    ptm <- proc.time() 
    fable_est <- FABLEPosteriorMean(data$Y)
    fable_time <- proc.time() - ptm
    fable_est$estRank
    fable_k_hat <- rbind(fable_k_hat, c(compute_metrics_fable(fable_est, compute_coverage=T, idx_cvg=idx_cvg, idx_fr=idx_fr),
                                        fable_time[3], fable_est$estRank))
    fable_k_hat
    k_jic <- fable_est$estRank
    rm(fable_est)
  }
  
  if(test_barigozzi_cho){
    # true k
    ptm <- proc.time() 
    bc_est <- barigozzi_cho_est(data$Y, k)
    bc_time <- proc.time() - ptm
    bc_true_k <- rbind(bc_true_k, c(compute_metrics(bc_est, idx_cvg=idx_cvg, idx_fr=idx_fr), bc_time[3], k))
    bc_true_k
    # k_hat b & n
    ptm <- proc.time() 
    bc_est <- barigozzi_cho_est(data$Y, k_bn$k_hat)
    bc_time <- proc.time() - ptm
    bc_k_bn <- rbind(bc_k_bn, c(compute_metrics(bc_est, idx_cvg=idx_cvg, idx_fr=idx_fr), bc_time[3], k_bn$k_hat))
    bc_k_bn
    # k_hat a & h
    ptm <- proc.time() 
    bc_est <- barigozzi_cho_est(data$Y, k_ah$k_hat)
    bc_time <- proc.time() - ptm
    bc_k_ah <- rbind(bc_k_ah, c(compute_metrics(bc_est, idx_cvg=idx_cvg, idx_fr=idx_fr), bc_time[3], k_ah$k_hat))
    # k_hat jic
    ptm <- proc.time() 
    bc_est <- barigozzi_cho_est(data$Y, k_jic)
    bc_time <- proc.time() - ptm
    bc_k_jic <- rbind(bc_k_jic, c(compute_metrics(bc_est, idx_cvg=idx_cvg, idx_fr=idx_fr), bc_time[3], k_jic))
    # k_hat jic + 5
    ptm <- proc.time() 
    bc_est <- barigozzi_cho_est(data$Y, k_jic+5)
    bc_time <- proc.time() - ptm
    bc_k_jic_over <- rbind(bc_k_jic_over, c(compute_metrics(bc_est, idx_cvg=idx_cvg, idx_fr=idx_fr), bc_time[3], k_jic +5))
    # k + 5
    ptm <- proc.time() 
    bc_est <- barigozzi_cho_est(data$Y, k_over)
    bc_time <- proc.time() - ptm
    bc_k_over <- rbind(bc_k_over, c(compute_metrics(bc_est, idx_cvg=idx_cvg, idx_fr=idx_fr), bc_time[3], k+5))
    rm(bc_est)
    
  }
  
  if(test_eigenbayes){
    ptm <- proc.time() 
    s_Y <- svd(data$Y)
    #U <- s_Y$u[,1:H]
    #V <- s_Y$v[,1:H]
    #d <- s_Y$d[1:H]
    #D <- diag(as.vector(d))
    #D_tilde <- diag(as.vector(d)) - mean(s_Y$d[-c(1:H)])
    # true k
    #ptm <- proc.time() 
    eb_est <- eigenbayes_point_est(data$Y, s_Y, k)
    eb_time <- proc.time() - ptm
    eb_true_k <- rbind(eb_true_k, c(compute_metrics(eb_est, compute_coverage=T, idx_cvg=idx_cvg, idx_fr=idx_fr), 
                                    eb_time[3], k))
    eb_true_k
    # k_hat b & n
    #ptm <- proc.time() 
    eb_est <- eigenbayes_point_est(data$Y, s_Y, k_bn$k_hat)
    #eb_time <- proc.time() - ptm
    eb_k_bn <- rbind(eb_k_bn, c(compute_metrics(eb_est, compute_coverage=T, idx_cvg=idx_cvg, idx_fr=idx_fr),
                                eb_time[3], k_bn$k_hat))
    # k_hat a & h
    #ptm <- proc.time() 
    eb_est <- eigenbayes_point_est(data$Y, s_Y, k_ah$k_hat)
    #eb_time <- proc.time() - ptm
    eb_k_ah <- rbind(eb_k_ah, c(compute_metrics(eb_est, compute_coverage=T, idx_cvg=idx_cvg, idx_fr=idx_fr),
                                eb_time[3], k_ah$k_hat))
    # k_hat jic
    #ptm <- proc.time() 
    eb_est <- eigenbayes_point_est(data$Y, s_Y, k_jic)
    #eb_time <- proc.time() - ptm
    eb_k_jic <- rbind(eb_k_jic, c(compute_metrics(eb_est, compute_coverage=T, idx_cvg=idx_cvg, idx_fr=idx_fr), 
                                  eb_time[3], k_jic))
    # k_hat jic + 5
    #ptm <- proc.time() 
    eb_est <- eigenbayes_point_est(data$Y, s_Y, k_jic+5)
    #eb_time <- proc.time() - ptm
    eb_k_jic_over <- rbind(eb_k_jic_over, c(compute_metrics(eb_est, compute_coverage=T, idx_cvg=idx_cvg, idx_fr=idx_fr),
                                            eb_time[3], k_jic +5))
    # k + 5
    #ptm <- proc.time() 
    eb_est <- eigenbayes_point_est(data$Y, s_Y, k_over)
    #eb_time <- proc.time() - ptm
    eb_k_over <- rbind(eb_k_over, c(compute_metrics(eb_est, compute_coverage=T, idx_cvg=idx_cvg, idx_fr=idx_fr),
                                    eb_time[3], k+5))
    rm(eb_est)
    
  }
  
  if(test_rotate){
    # true k
    ptm <- proc.time() 
    rot_est <- rotate_est(data$Y, k)
    rotate_time <- proc.time() - ptm
    rotate_true_k <- rbind(rotate_true_k, c(compute_metrics(rot_est, idx_cvg=idx_cvg, idx_fr=idx_fr), rotate_time[3], k))
    # k_hat b & n
    ptm <- proc.time() 
    rot_est <- rotate_est(data$Y, k_bn$k_hat)
    rotate_time <- proc.time() - ptm
    rotate_k_bn <- rbind(rotate_k_bn, c(compute_metrics(rot_est, idx_cvg=idx_cvg, idx_fr=idx_fr), rotate_time[3], k_bn$k_hat))
    # k_hat a & h
    ptm <- proc.time() 
    rot_est <- rotate_est(data$Y, k_ah$k_hat)
    rotate_time <- proc.time() - ptm
    rotate_k_ah <- rbind(rotate_k_ah, c(compute_metrics(rot_est, idx_cvg=idx_cvg, idx_fr=idx_fr), rotate_time[3], k_ah$k_hat))
    # k_hat jic
    ptm <- proc.time() 
    rot_est <- rotate_est(data$Y, k_jic)
    rot_time <- proc.time() - ptm
    rotate_k_jic <- rbind(rotate_k_jic, c(compute_metrics(rot_est, idx_cvg=idx_cvg, idx_fr=idx_fr), rot_time[3], k_jic))
    # k_hat jic + 5
    ptm <- proc.time() 
    rot_est <- rotate_est(data$Y, k_jic+5)
    rot_time <- proc.time() - ptm
    rotate_k_jic_over <- rbind(rotate_k_jic_over, c(compute_metrics(rot_est, idx_cvg=idx_cvg, idx_fr=idx_fr), rot_time[3], k_jic + 5))
    # k + 5
    ptm <- proc.time() 
    rot_est <- rotate_est(data$Y, k_over)
    rotate_time <- proc.time() - ptm
    rotate_k_over <- rbind(rotate_k_over, c(compute_metrics(rot_est,idx_cvg=idx_cvg, idx_fr=idx_fr), rotate_time[3], k+5))
    rm(rot_est)
  }
  
  if(test_pca){
    # true k
    ptm <- proc.time() 
    pca_est <- spectral_est(data$Y, s_Y, k)
    pca_time <- proc.time() - ptm
    pca_true_k <- rbind(pca_true_k, c(compute_metrics(pca_est, idx_cvg=idx_cvg, idx_fr=idx_fr), pca_time[3], k))
    # k_hat b & n
    ptm <- proc.time() 
    pca_est <- spectral_est(data$Y, s_Y, k_bn$k_hat)
    pca_time <- proc.time() - ptm
    pca_k_bn <- rbind(pca_k_bn, c(compute_metrics(pca_est, idx_cvg=idx_cvg, idx_fr=idx_fr), pca_time[3], k_bn$k_hat))
    # k_hat a & h. 
    ptm <- proc.time() 
    pca_est <- spectral_est(data$Y, s_Y, k_ah$k_hat)
    pca_time <- proc.time() - ptm
    pca_k_ah <- rbind(pca_k_ah, c(compute_metrics(pca_est, idx_cvg=idx_cvg, idx_fr=idx_fr), pca_time[3], k_ah$k_hat))
    # k_hat jic
    ptm <- proc.time() 
    pca_est <- spectral_est(data$Y, s_Y, k_jic)
    pca_time <- proc.time() - ptm
    pca_k_jic <- rbind(pca_k_jic, c(compute_metrics(pca_est, idx_cvg=idx_cvg, idx_fr=idx_fr), pca_time[3], k_jic))
    # k_hat jic + 5
    ptm <- proc.time() 
    pca_est <- spectral_est(data$Y, s_Y, k_jic+5)
    pca_time <- proc.time() - ptm
    pca_k_jic_over <- rbind(pca_k_jic_over, c(compute_metrics(pca_est, idx_cvg=idx_cvg, idx_fr=idx_fr), pca_time[3], k_jic + 5))
    # k + 5
    ptm <- proc.time() 
    pca_est <- spectral_est(data$Y, s_Y, k_over)
    pca_time <- proc.time() - ptm
    pca_k_over <- rbind(pca_k_over, c(compute_metrics(pca_est, idx_cvg=idx_cvg, idx_fr=idx_fr), pca_time[3], k_over))
  }
  
}

names <- c('L_fr', 'C_fr', 'time', k)


names(bc_true_k) <- names
names(bc_k_bn)<- names
names(bc_k_ah)<- names
names(bc_k_jic)<- names
names(bc_k_jic_over)<- names
names(bc_k_over) <- names

names(eb_true_k) <- names
names(eb_k_bn) <- names
names(eb_k_ah) <- names
names(eb_k_jic)<- names
names(eb_k_jic_over)<- names
names(eb_k_over)<- names

names(rotate_true_k) <- names
names(rotate_k_bn) <- names
names(rotate_k_ah) <- names
names(rotate_k_jic)<- names
names(rotate_k_jic_over)<- names
names(rotate_k_over) <- names

names(pca_true_k) <- names
names(pca_k_bn) <- names
names(pca_k_ah) <- names
names(pca_k_jic)<- names
names(pca_k_jic_over)<- names
names(pca_k_over) <- names

names(fable_k_hat) <- names


colMeans(bc_true_k)
colMeans(bc_k_bn)
colMeans(bc_k_ah)
colMeans(bc_k_jic)
colMeans(bc_k_jic_over)
colMeans(bc_k_over) 

colMeans(eb_true_k) 
colMeans(eb_k_bn) 
colMeans(eb_k_ah) 
colMeans(eb_k_jic)
colMeans(eb_k_jic_over)
colMeans(eb_k_over)

colMeans(rotate_true_k) 
colMeans(rotate_k_bn) 
colMeans(rotate_k_ah) 
colMeans(rotate_k_jic)
colMeans(rotate_k_jic_over)
colMeans(rotate_k_over) 

colMeans(pca_true_k) 
colMeans(pca_k_bn) 
colMeans(pca_k_ah) 
colMeans(pca_k_jic)
colMeans(pca_k_jic_over)
colMeans(pca_k_over) 

colMeans(fable_k_hat) 

set.seed(123)



write.csv(fama_results, paste0('simulations/results/scenario_', scenario, '/fama_results.csv'))
write.csv(fable_results, paste0('simulations/results/scenario_', scenario, '/fable_results.csv'))
write.csv(mofa_results, paste0('simulations/results/scenario_', scenario, '/mofa_results.csv'))
write.csv(rotate_results, paste0('simulations/results/scenario_', scenario, '/rotate_results.csv'))
save.image(file=paste0('simulations/results/scenario_', scenario, '/scenario_', scenario,'.RData'))

library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(patchwork)

fama_results <- fama_results %>% mutate(method = "FAMA")
fable_results <- fable_results %>% mutate(method = "FABLE")
mofa_results <- mofa_results %>% mutate(method = "MOFA")
rotate_results <- rotate_results %>% mutate(method = "ROTATE")

all_results <- bind_rows(fama_results, fable_results, mofa_results, rotate_results)
all_results$method <- factor(all_results$method, levels = c("FAMA", "FABLE", "MOFA", "ROTATE"))


p1 <- ggplot(all_results, aes(x = method, y = rmse_all, fill = method)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "RMSE Overall", y = "RMSE", x = "Method") +
  theme(legend.position = "none")

intra_data <- all_results %>%
  dplyr::select(method, starts_with("rmses_intra")) %>%
  pivot_longer(cols = starts_with("rmses_intra"), names_to = "metric", values_to = "value")

p2 <- ggplot(intra_data, aes(x = method, y = value, fill = method)) +
  geom_boxplot(position = position_dodge()) +
  theme_minimal() +
  labs(title = "RMSEs Intra-view", x = "Method", y = "RMSE") +
  theme(legend.position = "none")

inter_data <- all_results %>%
  dplyr::select(method, starts_with("rmses_inter")) %>%
  pivot_longer(cols = starts_with("rmses_inter"), names_to = "metric", values_to = "value")

p3 <- ggplot(inter_data, aes(x = method, y = value, fill = method)) +
  geom_boxplot(position = position_dodge()) +
  theme_minimal() +
  labs(title = "RMSEs Inter-view", x = "Method", y = "RMSE") +
  theme(legend.position = "none")

p4 <- ggplot(all_results, aes(x = method, y = time_pe, fill = method)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Running Time", x = "Method", y = "Seconds") +
  theme(legend.position = "none")


dev.new()
(p1 / p2 / p3 / p4) + plot_layout(ncol = 4)



rbind(colMeans(fama_results[ , grepl("len", names(fama_results))]),
      rep(colMeans(fable_results[ , grepl("len", names(fable_results))]),2))

colMeans(fama_results)
colMeans(fable_results)


coverage_fama <- fama_results[ , grepl("cov", names(fama_results))]
coverage_fama_clt <- coverage_fama[ , grepl("clt", names(coverage_fama))]
coverage_fama_bvm <- coverage_fama[ , grepl("posterior", names(coverage_fama))]

coverage_fable <- fable_results[ , grepl("cov", names(fable_results))]

coverage_fama_clt[] <- lapply(coverage_fama_clt, as.double)
coverage_fama_bvm[] <- lapply(coverage_fama_bvm, as.double)
coverage_fable[] <- lapply(coverage_fable, as.double)

# intra 
coverage_fama_clt_intra <- coverage_fama_clt[ , grepl("intra", names(coverage_fama_clt))]
coverage_fama_bvm_intra <- coverage_fama_bvm[ , grepl("intra", names(coverage_fama_bvm))]
coverage_fable_intra <- coverage_fable[ , grepl("intra", names(coverage_fable))]

mean(colMeans(coverage_fama_clt_intra));
mean(colMeans(coverage_fama_bvm_intra));
mean(colMeans(coverage_fable_intra));


# inter 
coverage_fama_clt_inter <- coverage_fama_clt[ , grepl("inter", names(coverage_fama_clt))]
coverage_fama_bvm_inter <- coverage_fama_bvm[ , grepl("inter", names(coverage_fama_bvm))]
coverage_fable_inter <- coverage_fable[ , grepl("inter", names(coverage_fable))]

mean(colMeans(coverage_fama_clt_inter));
mean(colMeans(coverage_fama_bvm_inter));
mean(colMeans(coverage_fable_inter))
