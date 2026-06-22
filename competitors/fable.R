
library(FABLE)
Rcpp::sourceCpp("competitors/FABLE/src/updated-FABLE-functions.cpp")

#sourceCpp('competitors/fable_helpers_updated.cpp')

# slight modification of code from https://github.com/shounakch/FABLE/tree/main
# to implement the FABLE methodology from
# Shounak Chattopadhyay, Anru R. Zhang, and David B. Dunson, (2024) 
# "Blessing of dimension in Bayesian inference on covariance matrices"
# arXiv precmprint arXiv:2404.03805




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



FABLEPosteriorMean_2 <- function(Y,
                                 gamma0 = 1,
                                 delta0sq = 1,
                                 maxProp = 0.95, k=NA) {
  
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
  
  if(is.na(k)){
    CovEstMod = CPPFABLEPostMean(Y, gamma0, delta0sq, U_Y, V_Y, svalsY, kMax)
  } else {
    CovEstMod = CPPFABLEPostMean_updated(Y, gamma0, delta0sq, U_Y, V_Y, svalsY, kMax, kInput=k )
  }
  CovEst = CovEstMod$CovPostMean
  kEst = CovEstMod$kEst
  
  FABLEHypPars = FABLEHyperParameters(Y, U_Y, V_Y, svalsY, kEst)
  
  Part1 = FABLEHypPars$G
  #print(dim(Part1))
  # Part2 = as.numeric(FABLEHypPars$gammaDeltasq / (FABLEHypPars$gamman - 2))
  # CovEst = Part1 + diag(Part2)
  
  tFABLEPostMean2 = proc.time()
  tPostMean = (tFABLEPostMean2 - tFABLEPostMean1)[3]
  
  OutputList = list("FABLEPostMean" = CovEst,
                    "Lambda_outer" = Part1,
                    "svdY" = svdY,
                    "estRank" = kEst,
                    "runTime" = tPostMean)
  
  return(OutputList)
  
}


FABLEPosteriorSampler_2 <- function(Y,
                                  gamma0 = 1,
                                  delta0sq = 1,
                                  maxProp = 0.95,
                                  MC = 1000, k=NA) {
  
  tFABLESample1 = proc.time()
  
  Y = as.matrix(Y)
  n = nrow(Y)
  p = ncol(Y)
  svdY = svd(Y)
  U_Y = svdY$u
  V_Y = svdY$v
  svalsY = svdY$d
  kMax = min(which(cumsum(svalsY) / sum(svalsY) >= maxProp))
  
  if(is.na(k)){
    kEst = CPPRankEstimator(Y, 
                            U_Y,
                            V_Y,
                            svalsY,
                            kMax)
  } else {
    kEst = k
  }
  
  
  FABLEHypPars = FABLEHyperParameters(Y,
                                      U_Y,
                                      V_Y,
                                      svalsY,
                                      kEst)
  
  CovCorrectMatrix = cov_correct_matrix(FABLEHypPars$SigmaSqEstimate, 
                                        FABLEHypPars$G)
  
  varInflation = (sum(CovCorrectMatrix) / (p*(p+1)/2))^2
  
  FABLESamples = CPPFABLESampler(Y, 
                                 gamma0, 
                                 delta0sq, 
                                 MC,
                                 U_Y,
                                 V_Y,
                                 svalsY,
                                 kEst,
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

