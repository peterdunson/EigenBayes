ks <- seq(10, 150, by=10)

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

fable_log_lik <- dmvnorm(Y_test_gene,  mu=rep(0, p), Sigma=fable_est$FABLEPostMean, log=T)
sum(fable_log_lik)


for(k in ks){
  
  print(k)
  
  eb_ests[[as.character(k)]] <- eigenbayes_point_est(Y_train_gene, s_Y, k)
  rot_ests[[as.character(k)]] <- rotate_est(Y_train_gene, k)
  pca_ests[[as.character(k)]] <- spectral_est(Y_train_gene, s_Y, k)
  bc_ests[[as.character(k)]] <- barigozzi_cho_est(Y_train_gene, k)
  
  eb_ests[[as.character(k)]]$cov_hat <- tcrossprod(eb_ests[[as.character(k)]]$Lambda_hat) + diag(eb_ests[[as.character(k)]]$sigmas_sq_mean)
  eb_log_liks[[as.character(k)]] <- dmvnorm(Y_test_gene,  mu=rep(0, p), Sigma=eb_ests[[as.character(k)]]$cov_hat, log=T)
  print(sum(eb_log_liks[[as.character(k)]]))
  
  rot_ests[[as.character(k)]]$cov_hat <- tcrossprod(rot_ests[[as.character(k)]]$Lambda_hat) + diag(rot_ests[[as.character(k)]]$sigmas_sq_hat)
  rot_log_liks[[as.character(k)]] <- dmvnorm(Y_test_gene,  mu=rep(0, p), Sigma=rot_ests[[as.character(k)]]$cov_hat, log=T)
  print(sum(rot_log_liks[[as.character(k)]]))
  
  pca_ests[[as.character(k)]]$cov_hat <- tcrossprod(pca_ests[[as.character(k)]]$Lambda_hat) + diag(pca_ests[[as.character(k)]]$sigmas_sq_hat)
  pca_log_liks[[as.character(k)]] <- dmvnorm(Y_test_gene,  mu=rep(0, p), Sigma=pca_ests[[as.character(k)]]$cov_hat, log=T)
  print(sum(pca_log_liks[[as.character(k)]]))
  
  bc_ests[[as.character(k)]]$cov_hat <- tcrossprod(bc_ests[[as.character(k)]]$Lambda_hat) + diag(bc_ests[[as.character(k)]]$sigmas_sq_hat)
  bc_log_liks[[as.character(k)]] <- dmvnorm(Y_test_gene,  mu=rep(0, p), Sigma=bc_ests[[as.character(k)]]$cov_hat, log=T)
  print(sum(bc_log_liks[[as.character(k)]]))
  
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


pval_res <- data.frame(
  k = ks,
  EB_vs_BC = NA,
  EB_vs_PCA = NA,
  EB_vs_ROTATE = NA,
  EB_vs_FABLE = NA
)

for(k in ks){
  
  test.1 <- t.test(eb_log_liks[[as.character(k)]], bc_log_liks[[as.character(k)]], alternative='greater', paired=TRUE)
  pval_res[pval_res$k == k, "EB_vs_BC"] <- test.1$p.value
  
  test.1 <- t.test(eb_log_liks[[as.character(k)]], pca_log_liks[[as.character(k)]], alternative='greater', paired=TRUE)
  pval_res[pval_res$k == k, "EB_vs_PCA"] <- test.1$p.value
  
  test.1 <- t.test(eb_log_liks[[as.character(k)]], rot_log_liks[[as.character(k)]], alternative='greater', paired=TRUE)
  pval_res[pval_res$k == k, "EB_vs_ROTATE"] <- test.1$p.value
  
  test.1 <- t.test(eb_log_liks[[as.character(k)]], fable_log_lik, alternative='greater', paired=TRUE)
  pval_res[pval_res$k == k, "EB_vs_FABLE"] <- test.1$p.value
  
}

pval_res

library(ggplot2)
library(tidyr)

oos_log_lik_plot <- pivot_longer(
  oos_log_lik_res,
  cols = c("EB", "BC", "PCA", "ROTATE", "FABLE"),
  names_to = "method",
  values_to = "oos_log_lik"
)

ggplot(oos_log_lik_plot, aes(x=k, y=oos_log_lik, group=method, linetype=method, shape=method)) +
  geom_line() +
  geom_point(size=2) +
  theme_bw() +
  labs(x="k", y="OOS log-likelihood")