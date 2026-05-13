
library(FABLE)

# slight modification of code from https://github.com/shounakch/FABLE/tree/main
# to implement the FABLE methodology from
# Shounak Chattopadhyay, Anru R. Zhang, and David B. Dunson, (2024) 
# "Blessing of dimension in Bayesian inference on covariance matrices"
# arXiv precmprint arXiv:2404.03805


construct_fable_cov_samples <- function(FABLESamples, idx=NULL){
  n_MC <- nrow(FABLESamples$CCFABLESamples$LambdaSamples)
  k <- FABLESamples$estRank
  
  
  p_long <- nrow(matrix(FABLESamples$CCFABLESamples$LambdaSamples[1,], ncol = k, byrow=T))
  lr_samples <- array(NA, dim=c(n_MC, p_long, p_long))
  
  p <- length(idx)
  cov_samples <- array(NA, dim=c(n_MC, p, p))
  
  for(t in 1:n_MC){
    Lambda_t <- (matrix(FABLESamples$CCFABLESamples$LambdaSamples[t,], ncol = k, byrow=T))
    lr_samples[t,,] <- tcrossprod(Lambda_t)
    
    cov_samples[t,,] <- lr_samples[t,idx,idx]  +
      diag(FABLESamples$CCFABLESamples$SigmaSqSamples[t,idx])
  }
  cov_mean <- apply(cov_samples, c(2,3), mean)
  Lambda_outer_mean <- apply(lr_samples, c(2,3), mean)
  
  return(list(cov_mean=cov_mean, cov_samples=cov_samples, Lambda_outer_mean=Lambda_outer_mean))
  
}


fable_low_rank_signal <- function(Y, FABLESamples){
  n_MC <- nrow(FABLESamples$CCFABLESamples$LambdaSamples)
  k <- FABLESamples$estRank
  p <- ncol(Y)
  n <- nrow(Y)
  
  lr_signal <- matrix(0, n, p)
  
  for(t in 1:n_MC){
    Lambda_t <- (matrix(FABLESamples$CCFABLESamples$LambdaSamples[t,], ncol = k, byrow=T))
    Sigma_t <- FABLESamples$CCFABLESamples$SigmaSqSamples[t,]
    var_t <- solve(t(Lambda_t) %*% diag(1/Sigma_t) %*% Lambda_t)
    factors_mean <- latent_factors_full_conditional_mean(
      Y, Lambda_t, Sigma_t, var_t
    )
    lr_signal <- lr_signal + factors_mean %*% t(Lambda_t)
  }
  lr_signal <- lr_signal / n_MC
  return(lr_signal)
  
}






fable_monte_carlo_ci <- function(fable_cov_samples, alpha=0.05){
  idx <- upper.tri(fable_cov_samples$cov_mean, diag=T)
  cis <- apply(fable_cov_samples$cov_samples, c(2,3),function(x)(quantile(x, probs=c(alpha/2, 1-alpha/2))))
  return(list(credible_intervals=cis))
}


PseudoPosteriorMean_3 <- function(Y,
                                gamma0 = 1,
                                delta0sq = 1,
                                maxProp = 0.5) {
  
  tFABLEPostMean1 = proc.time()
  
  Y = as.matrix(Y)
  n = nrow(Y)
  p = ncol(Y)
  svdY = svd(Y)
  U_Y = svdY$u
  V_Y = svdY$v
  svalsY = svdY$d
  kMax = min(which(cumsum(svalsY) / sum(svalsY) >= maxProp))
  
  kEst = RankEstimator(Y, 
                       U_Y,
                       V_Y,
                       svalsY,
                       kMax)
  
  FABLEHypPars = FABLEHyperParameters(Y,
                                      U_Y,
                                      V_Y,
                                      svalsY,
                                      kEst,
                                      gamma0,
                                      delta0sq)
  
  Part1 = FABLEHypPars$G
  Part2 = as.numeric(FABLEHypPars$gammaDeltasq / (FABLEHypPars$gamman - 2))
  CovEst = Part1 + diag(Part2)
  
  tFABLEPostMean2 = proc.time()
  tPostMean = (tFABLEPostMean2 - tFABLEPostMean1)[3]
  
  
  M_hat = U_Y[,1:kEst, drop=FALSE] * sqrt(n)
  Lambda_hat = V_Y[,1:kEst, drop=FALSE] %*% diag(svdY$d[1:kEst, drop=F]) * sqrt(n) / (n + 1/FABLEHypPars$tausq_est)
  print(dim(M_hat))
  print(dim(Lambda_hat))
  print(FABLEHypPars$tausq_est)
  OutputList = list("FABLEPostMean" = CovEst,
                    "Lambda_outer" = Part1,
                    "FABLEHyperParameters" = FABLEHypPars,
                    "svdY" = svdY,
                    "estRank" = kEst,
                    "runTime" = tPostMean,
                    'X' = M_hat %*% t(Lambda_hat))
  
  return(OutputList)
  
}



FABLEPosteriorMean_2 <- function(Y,
                               gamma0 = 1,
                               delta0sq = 1,
                               maxProp = 0.95) {
  
  tFABLEPostMean1 = proc.time()
  
  Y = as.matrix(Y)
  n = nrow(Y)
  p = ncol(Y)
  svdY = svd(Y)
  U_Y = svdY$u
  V_Y = svdY$v
  svalsY = svdY$d
  kMax = min(which(cumsum(svalsY) / sum(svalsY) >= maxProp))
  
  # kEst = RankEstimator(Y, 
  #                      U_Y,
  #                      V_Y,
  #                      svalsY,
  #                      kMax)
  
  CovEstMod = CPPFABLEPostMean(Y, gamma0, delta0sq, U_Y, V_Y, svalsY, kMax)
  CovEst = CovEstMod$CovPostMean
  kEst = CovEstMod$kEst
  
  # FABLEHypPars = FABLEHyperParameters(Y,
  #                                     U_Y,
  #                                     V_Y,
  #                                     svalsY,
  #                                     kEst,
  #                                     gamma0,
  #                                     delta0sq)
  
  # Part1 = FABLEHypPars$G
  # Part2 = as.numeric(FABLEHypPars$gammaDeltasq / (FABLEHypPars$gamman - 2))
  # CovEst = Part1 + diag(Part2)
  
  tFABLEPostMean2 = proc.time()
  tPostMean = (tFABLEPostMean2 - tFABLEPostMean1)[3]
  
  
  #sigsq_hat_diag = (sum(square(Y - UDVt), 0)).t() / n; 
  sigsq_hat_diag = colMeans((Y - tcrossprod(U_Y[,1:kEst, drop=FALSE])%*%Y)^2)
  #tausq_est = (mean(sum(square(YtU.t()), 0).t() / sigsq_hat_diag)) / (n * k);
  tausq_est = (mean(rowSums((crossprod(Y, U_Y[,1:kEst, drop=FALSE]))^2) / sigsq_hat_diag)) / (n * kEst);
  M_hat = U_Y[,1:kEst, drop=FALSE] * sqrt(n)
  Lambda_hat = V_Y[,1:kEst, drop=FALSE] %*% diag(svdY$d[1:kEst, drop=F]) * sqrt(n) / (n + 1/tausq_est)
  print(tausq_est)
  
  
  OutputList = list("FABLEPostMean" = CovEst,
                    "Lambda_outer" = tcrossprod(Lambda_hat),
                    "svdY" = svdY,
                    "estRank" = kEst,
                    "runTime" = tPostMean,
                    'Y_hat' = M_hat %*% t(Lambda_hat))
  
  return(OutputList)
  
}



PseudoPosteriorSampler_2 <- function(fit,
                                   Y,
                                   gamma0 = 1,
                                   delta0sq = 1,
                                   maxProp = 0.5,
                                   MC = 1000) {
  
  tFABLESample1 = proc.time()
  
  Y = as.matrix(Y)
  n = nrow(Y)
  p = ncol(Y)
  svdY = svd(Y)
  U_Y = svdY$u
  V_Y = svdY$v
  svalsY = svdY$d
  #kMax = min(which(cumsum(svalsY) / sum(svalsY) >= maxProp))
  kEst = fit$estRank
  
  FABLEHypPars = FABLEHyperParameters(Y,
                                      U_Y,
                                      V_Y,
                                      svalsY,
                                      kEst,
                                      gamma0,
                                      delta0sq)
  
  CovCorrectMatrix = cov_correct_matrix(FABLEHypPars$SigmaSqEstimate, 
                                        FABLEHypPars$G)
  
  varInflation = (sum(CovCorrectMatrix) / (p*(p+1)/2))^2
  
  FABLESamples = FABLESampler(Y, 
                              gamma0, 
                              delta0sq, 
                              MC,
                              U_Y,
                              V_Y,
                              svalsY,
                              kEst,
                              FABLEHypPars$tauSqEstimate,
                              FABLEHypPars$gammaDeltasq,
                              FABLEHypPars$G0,
                              varInflation)
  
  tFABLESample2 = proc.time()
  tSample = (tFABLESample2 - tFABLESample1)[3]
  
  OutputList = list("CCFABLESamples" = FABLESamples,
                    "FABLEHyperParameters" = FABLEHypPars,
                    "svdY" = svdY,
                    "estRank" = kEst,
                    "varInflation" = varInflation,
                    "runTime" = tSample)
  
  return(OutputList)
  
}

