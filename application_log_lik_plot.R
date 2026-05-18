


source('eigenbayes_functions.R')


Y_norm <- readRDS('data/Y_gene_norm.rds')
n <- nrow(Y_norm)

id = as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))
print(id)
set.seed(id)
test_gene <- sample(1:n, floor(0.2 * n))
Y_train_gene <- as.matrix(Y_norm[-test_gene, ])
Y_test_gene <- as.matrix(Y_norm[test_gene, ])



ks <- seq(40, 160, by=20)

s_Y <- svd(Y_train_gene)

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


library(emdbook)

p <- ncol(Y_test_gene)
fable_est <- FABLEPosteriorMean(Y_train_gene)
k_jic <- fable_est$estRank
k_over <- floor(sqrt(n))
k_over_2 <- floor(sqrt(p))
fable_log_lik <- dmvnorm(Y_test_gene,  mu=rep(0, p), Sigma=fable_est$FABLEPostMean, log=T)
sum(fable_log_lik)



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
  FABLE = rep(sum(fable_log_lik), length(ks))
)

oos_log_lik_res

saveRDS(oos_log_lik_res, file = paste0("results/application/", id, ".rds"))


