


source('application/preprocess_application.R')



n <- nrow(Y_norm)

# iterarate for id in 1:50
id = 1
print(id)
set.seed(id)
test_gene <- sample(1:n, floor(0.2 * n))
Y_train_gene <- as.matrix(Y_norm[-test_gene, ])
Y_test_gene <- as.matrix(Y_norm[test_gene, ])



ks <- seq(40, 160, by=10)


eb_ests <- vector("list", length(ks))
rot_ests <- vector("list", length(ks))
pca_ests <- vector("list", length(ks))
bc_ests <- vector("list", length(ks))

names(eb_ests) <- ks
names(rot_ests) <- ks
names(pca_ests) <- ks
names(bc_ests) <- ks

eb_log_liks <- vector("list", length(ks))
rot_log_liks <- vector("list", length(ks))
pca_log_liks <- vector("list", length(ks))
bc_log_liks <- vector("list", length(ks))

names(eb_log_liks) <- ks
names(rot_log_liks) <- ks
names(pca_log_liks) <- ks
names(bc_log_liks) <- ks



s_Y <- svd(Y_train_gene)


p <- ncol(Y_test_gene)
fable_est <- FABLEPosteriorMean(Y_train_gene)
k_jic <- fable_est$estRank
k_over <- floor(sqrt(n))
k_over_2 <- floor(sqrt(p))
fable_log_lik <- dmvnorm(Y_test_gene,  mu=rep(0, p), Sigma=fable_est$FABLEPostMean, log=T)
sum(fable_log_lik)

ptr <- proc.time()
out_msgp = linearMGSP(X = Y_train_gene, nrun = 3000, burn = 1000, output=c('covMean', 'numFactors'), verbose=T)
mgsp_time <- proc.time() - ptr
mgsp_time[3]



for(k in ks){
  
  print(k)
  
  eb_est <- eigenbayes_point_est(Y_train_gene, s_Y, k)
  eb_cov_hat <- tcrossprod(eb_est$Lambda_hat) + diag(eb_est$sigmas_sq_mean)
  eb_log_liks[[as.character(k)]] <- dmvnorm(Y_test_gene,  mu=rep(0, p), Sigma=eb_cov_hat, log=T)
  print(sum(eb_log_liks[[as.character(k)]]))
  rm(eb_est, eb_cov_hat)
  gc()
  
  rot_est <- rotate_est(Y_train_gene, k)
  rot_cov_hat <- tcrossprod(rot_est$Lambda_hat) + diag(rot_est$sigmas_sq_hat)
  rot_log_liks[[as.character(k)]] <- dmvnorm(Y_test_gene,  mu=rep(0, p), Sigma=rot_cov_hat, log=T)
  print(sum(rot_log_liks[[as.character(k)]]))
  rm(rot_est, rot_cov_hat)
  gc()
  
  pca_est <- spectral_est(Y_train_gene, s_Y, k)
  pca_cov_hat <- tcrossprod(pca_est$Lambda_hat) + diag(pca_est$sigmas_sq_hat)
  pca_log_liks[[as.character(k)]] <- dmvnorm(Y_test_gene,  mu=rep(0, p), Sigma=pca_cov_hat, log=T)
  print(sum(pca_log_liks[[as.character(k)]]))
  rm(pca_est, pca_cov_hat)
  gc()
  
  bc_est <- barigozzi_cho_est(Y_train_gene, k)
  bc_cov_hat <- tcrossprod(bc_est$Lambda_hat) + diag(bc_est$sigmas_sq_hat)
  bc_log_liks[[as.character(k)]] <- dmvnorm(Y_test_gene,  mu=rep(0, p), Sigma=bc_cov_hat, log=T)
  print(sum(bc_log_liks[[as.character(k)]]))
  rm(bc_est, bc_cov_hat)
  gc()
  
}

oos_log_lik_res <- data.frame(
  k = ks,
  EB = sapply(eb_log_liks, sum),
  BC = sapply(bc_log_liks, sum),
  PCA = sapply(pca_log_liks, sum),
  ROTATE = sapply(rot_log_liks, sum),
  FABLE = rep(sum(fable_log_lik), length(ks)),
  MGSP = rep(sum(mgsp_log_lik), length(ks))
)

oos_log_lik_res
#saveRDS(oos_log_lik_res, file = paste0("results/application/", id, ".rds"))

extra_ks <- c(k_jic + 10, k_over_2)
names(extra_ks) <- c("JIC+10", "sqrt(p)")
extra_log_liks <- list()

for(k_name in names(extra_ks)){
  
  k <- extra_ks[[k_name]]
  print(k_name)
  print(k)
  
  eb_est <- eigenbayes_point_est(Y_train_gene, s_Y, k)
  eb_cov_hat <- tcrossprod(eb_est$Lambda_hat) + diag(eb_est$sigmas_sq_mean)
  eb_log_lik <- dmvnorm(Y_test_gene,  mu=rep(0, p), Sigma=eb_cov_hat, log=T)
  print(sum(eb_log_lik))
  rm(eb_est, eb_cov_hat)
  gc()
  
  rot_est <- rotate_est(Y_train_gene, k)
  rot_cov_hat <- tcrossprod(rot_est$Lambda_hat) + diag(rot_est$sigmas_sq_hat)
  rot_log_lik <- dmvnorm(Y_test_gene,  mu=rep(0, p), Sigma=rot_cov_hat, log=T)
  print(sum(rot_log_lik))
  rm(rot_est, rot_cov_hat)
  gc()
  
  pca_est <- spectral_est(Y_train_gene, s_Y, k)
  pca_cov_hat <- tcrossprod(pca_est$Lambda_hat) + diag(pca_est$sigmas_sq_hat)
  pca_log_lik <- dmvnorm(Y_test_gene,  mu=rep(0, p), Sigma=pca_cov_hat, log=T)
  print(sum(pca_log_lik))
  rm(pca_est, pca_cov_hat)
  gc()
  
  bc_est <- barigozzi_cho_est(Y_train_gene, k)
  bc_cov_hat <- tcrossprod(bc_est$Lambda_hat) + diag(bc_est$sigmas_sq_hat)
  bc_log_lik <- dmvnorm(Y_test_gene,  mu=rep(0, p), Sigma=bc_cov_hat, log=T)
  print(sum(bc_log_lik))
  rm(bc_est, bc_cov_hat)
  gc()
  
  extra_log_liks[[k_name]] <- list(
    k = k,
    EB = eb_log_lik,
    BC = bc_log_lik,
    PCA = pca_log_lik,
    ROTATE = rot_log_lik,
    FABLE = fable_log_lik,
    MGSP = mgsp_log_lik
  )
  
  rm(eb_log_lik, bc_log_lik, pca_log_lik, rot_log_lik)
  gc()
  
}

extra_oos_log_lik_res <- do.call(rbind, lapply(names(extra_log_liks), function(k_name){
  
  data.frame(
    k_type = k_name,
    k = extra_log_liks[[k_name]]$k,
    EB = sum(extra_log_liks[[k_name]]$EB),
    BC = sum(extra_log_liks[[k_name]]$BC),
    PCA = sum(extra_log_liks[[k_name]]$PCA),
    ROTATE = sum(extra_log_liks[[k_name]]$ROTATE),
    FABLE = sum(extra_log_liks[[k_name]]$FABLE),
    MGSP = sum(extra_log_liks[[k_name]]$MGSP)
  )
  
}))

extra_pvals_res <- do.call(rbind, lapply(names(extra_log_liks), function(k_name){
  
  data.frame(
    k_type = k_name,
    k = extra_log_liks[[k_name]]$k,
    comparison = c("EB > BC", "EB > PCA", "EB > ROTATE", "EB > FABLE", "EB > MGSP" ),
    p_value = c(
      t.test(extra_log_liks[[k_name]]$EB, extra_log_liks[[k_name]]$BC, paired=TRUE, alternative="greater")$p.value,
      t.test(extra_log_liks[[k_name]]$EB, extra_log_liks[[k_name]]$PCA, paired=TRUE, alternative="greater")$p.value,
      t.test(extra_log_liks[[k_name]]$EB, extra_log_liks[[k_name]]$ROTATE, paired=TRUE, alternative="greater")$p.value,
      t.test(extra_log_liks[[k_name]]$EB, extra_log_liks[[k_name]]$FABLE, paired=TRUE, alternative="greater")$p.value,
      t.test(extra_log_liks[[k_name]]$EB, extra_log_liks[[k_name]]$MGSP, paired=TRUE, alternative="greater")$p.value
      
    )
  )
  
}))

extra_oos_log_lik_res
extra_pvals_res


#saveRDS(
#  extra_oos_log_lik_res,
#  file = paste0("results/application/extra_oos_log_lik_res_", id, ".rds")
#)

#saveRDS(
#  extra_pvals_res,
#  file = paste0("results/application/extra_pvals_res_", id, ".rds")
#)


rm(list=ls())


