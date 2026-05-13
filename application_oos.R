#genedata=readRDS('../BLAST/application/data/genedata.rds')

download_genedata=function(){
  library(httr)
  url <- "https://utdallas.box.com/shared/static/tuqwc8i0mzixs83wkvtla7qx0sg34365.rda"
  temp <- tempfile(fileext = ".rda")
  
  # Download the file to a temporary location
  httr::GET(url, write_disk(temp, overwrite = TRUE))
  
  # Load the file from the temporary location
  load(temp)
  
  # Clean up the temporary file
  unlink(temp)
  
  return(genedata)
}



genedata=download_genedata()

set.seed(123)

library(Biobase)
library(genefilter)
library(ggplot2)



array.data= genedata$array2
Ident.array= genedata$array2.types

cutoff <- 0.50
dat=ExpressionSet(assayData = array.data)
filter.dat=varFilter(dat, var.cutoff = cutoff)
hvgs_in_array=rownames(filter.dat)
rm(filter.dat,dat)

genes.use=hvgs_in_array
length(genes.use)

Y=t(array.data[genes.use,])

genes.common <- hvgs_in_array
set.seed(123)  # for reproducibility
genes.use <- sample(genes.common, 5000)

length(genes.use)

Y=t(array.data[genes.use,])

inv_cdf_normalize_same_sd <- function(x) {
  s <- sd(x, na.rm = TRUE)
  m <- mean(x, na.rm = TRUE)
  
  r <- rank(x, na.last = "keep", ties.method = "average")
  n <- sum(!is.na(x))
  z <- qnorm((r - 0.5) / n)
  
  z <- z - mean(z, na.rm = TRUE)
  z <- z / sd(z, na.rm = TRUE) * s
  
  z
}

invnorm_fit <- function(x_train) {
  x_train_obs <- x_train[!is.na(x_train)]
  n <- length(x_train_obs)
  
  list(
    train_sorted = sort(x_train_obs),
    mean_train   = mean(x_train_obs),
    sd_train     = sd(x_train_obs),
    n_train      = n
  )
}

invnorm_apply <- function(x, fit) {
  s <- fit$sd_train
  train_sorted <- fit$train_sorted
  n <- fit$n_train
  
  # empirical CDF based on training sample only
  p <- findInterval(x, train_sorted, left.open = FALSE) / n
  
  # avoid exactly 0 or 1
  p <- pmax(pmin(p, 1 - 1/(2*n)), 1/(2*n))
  
  z <- qnorm(p)
  
  # standardize using current mapped sample, then rescale to training sd
  z <- z - mean(z, na.rm = TRUE)
  z_sd <- sd(z, na.rm = TRUE)
  if (!is.finite(z_sd) || z_sd == 0) {
    z <- rep(0, length(z))
  } else {
    z <- z / z_sd * s
  }
  
  z
}

fit_apply_matrix_invnorm <- function(Y_train, Y_test) {
  p <- ncol(Y_train)
  
  Y_train_new <- matrix(NA_real_, nrow(Y_train), p)
  Y_test_new  <- matrix(NA_real_, nrow(Y_test),  p)
  
  colnames(Y_train_new) <- colnames(Y_train)
  colnames(Y_test_new)  <- colnames(Y_test)
  rownames(Y_train_new) <- rownames(Y_train)
  rownames(Y_test_new)  <- rownames(Y_test)
  
  fits <- vector("list", p)
  
  for (j in 1:p) {
    fits[[j]] <- invnorm_fit(Y_train[, j])
    Y_train_new[, j] <- invnorm_apply(Y_train[, j], fits[[j]])
    Y_test_new[, j]  <- invnorm_apply(Y_test[, j],  fits[[j]])
  }
  
  list(
    Y_train = Y_train_new,
    Y_test  = Y_test_new,
    fits    = fits
  )
}

n <- nrow(Y)
set.seed(123)
test_gene <- sample(1:n, floor(0.2 * n))

Y_train_gene <- as.matrix(Y[-test_gene, ])
Y_test_gene <- as.matrix(Y[test_gene, ])

tr1 <- fit_apply_matrix_invnorm(Y_train_gene, Y_test_gene)


gaussianize_train_test <- function(Y_train, Y_test, eps = 1e-6) {
  n <- nrow(Y_train)
  p <- ncol(Y_train)
  
  Y_train_g <- matrix(NA, n, p)
  Y_test_g  <- matrix(NA, nrow(Y_test), p)
  
  for (j in 1:p) {
    x_train <- Y_train[, j]
    
    # empirical CDF via ranks
    ranks <- rank(x_train, ties.method = "average")
    u_train <- (ranks - 0.5) / n
    
    # avoid infinities
    u_train <- pmin(pmax(u_train, eps), 1 - eps)
    
    # Gaussianized train
    Y_train_g[, j] <- qnorm(u_train)
    
    # store sorted values for interpolation
    x_sorted <- sort(x_train)
    u_sorted <- (1:n - 0.5) / n
    
    # map test using same CDF (interpolation)
    u_test <- approx(
      x = x_sorted,
      y = u_sorted,
      xout = Y_test[, j],
      rule = 2  # extrapolate at boundaries
    )$y
    
    u_test <- pmin(pmax(u_test, eps), 1 - eps)
    
    Y_test_g[, j] <- qnorm(u_test)
  }
  
  colnames(Y_train_g) <- colnames(Y_train)
  colnames(Y_test_g)  <- colnames(Y_test)
  
  list(train = Y_train_g, test = Y_test_g)
}



res <- gaussianize_train_test(Y_train_gene, Y_test_gene)
Y_train_gene <- res$train
Y_test_gene  <- res$test

hist(Y_train_gene[, 2])
hist(Y_test_gene[, 2])



p <- ncol(Y_test_gene)
set.seed(123)
impute_factor_index_gene <- sample(1:p, floor(0.5*p))
p_test <- p-length(impute_factor_index_gene)

fable_est <- FABLEPosteriorMean(Y_train_gene)
k_jic <- fable_est$estRank
k_over <- floor(sqrt(n))
k_over_2 <- floor(sqrt(p))

s_Y <- svd(Y_train_gene)
eb_est_k_jic <- eigenbayes_point_est(Y_train_gene, s_Y, k_jic)
eb_est_k_jic_over <- eigenbayes_point_est(Y_train_gene, s_Y, k_jic+10)
#eb_est_k_over <- eigenbayes_point_est(Y_train_gene, s_Y, k_over)
#eb_est_k_over_2 <- eigenbayes_point_est(Y_train_gene, s_Y, k_over_2)

rot_est_k_jic <- rotate_est(Y_train_gene, k_jic)
rot_est_k_jic_over <- rotate_est(Y_train_gene, k_jic+10)
#rot_est_k_over <- rotate_est(Y_train_gene, k_over)
#rot_est_k_over_2 <- rotate_est(Y_train_gene, k_over_2)

pca_est_k_jic <- spectral_est(Y_train_gene, s_Y, k_jic)
pca_est_k_jic_over <- spectral_est(Y_train_gene, s_Y, k_jic+5)

bc_est_k_jic <- barigozzi_cho_est(Y_train_gene, k_jic)
bc_est_k_jic_over <- barigozzi_cho_est(Y_train_gene, k_jic+5)


library(emdbook)

fable_log_lik <- dmvnorm(Y_test_gene,  mu=rep(0, p), Sigma=fable_est$FABLEPostMean, log=T)
sum(fable_log_lik)

eb_est_k_jic$cov_hat <- tcrossprod(eb_est_k_jic$Lambda_hat) + diag(eb_est_k_jic$sigmas_sq_mean)
eb_k_jic_log_lik <- dmvnorm(Y_test_gene,  mu=rep(0, p), Sigma=eb_est_k_jic$cov_hat, log=T)
sum(eb_k_jic_log_lik)
eb_est_k_jic_over$cov_hat <- tcrossprod(eb_est_k_jic_over$Lambda_hat) + diag(eb_est_k_jic_over$sigmas_sq_mean)
eb_k_jic_over_log_lik <- dmvnorm(Y_test_gene,  mu=rep(0, p), Sigma=eb_est_k_jic_over$cov_hat, log=T)

sum(eb_k_jic_over_log_lik)
sum(bc_k_jic_over_log_lik)
sum(pca_k_jic_over_log_lik)
sum(rot_k_jic_over_log_lik)
sum(fable_log_lik)


test.1 <- t.test(eb_k_jic_over_log_lik, bc_k_jic_over_log_lik, alternative='greater', paired=TRUE)
print(test.1$p.value)
test.1 <- t.test(eb_k_jic_over_log_lik, pca_k_jic_over_log_lik, alternative='greater', paired=TRUE)
print(test.1$p.value)
test.1 <- t.test(eb_k_jic_over_log_lik, rot_k_jic_over_log_lik, alternative='greater', paired=TRUE)
print(test.1$p.value)
test.1 <- t.test(eb_k_jic_over_log_lik, fable_log_lik, alternative='greater', paired=TRUE)
print(test.1$p.value)



test.1 <- t.test(eb_k_jic_over_log_lik, bc_k_jic_over_log_lik, alternative='greater', paired=TRUE)
print(test.1$p.value)
test.1 <- t.test(eb_k_jic_over_log_lik, bc_k_jic_over_log_lik, alternative='greater', paired=TRUE)
print(test.1$p.value)


predict_factors <- function(Y_old, Lambda_old, Psi_old, coverage){
  k <- ncol(Lambda_old)
  Var_fact <- solve(diag(1, k, k) + t(Lambda_old) %*% diag(1/Psi_old) %*% Lambda_old)
  mean_fact <- Y_old %*% diag(1/Psi_old) %*% Lambda_old  %*% Var_fact 
  if(coverage){  eta <-mean_fact + mvtnorm::rmvnorm(nrow(Y_old), sigma=Var_fact)}
  output <- list(Etas_mean=mean_fact)
  if(coverage){output$Etas_sample=eta}
  return(output)
}

predict_y <- function(Y_old, Lambda_new, Lambda_old, Psi_new, Psi_old, coverage=F){
  Eta <- predict_factors(Y_old, Lambda_old, Psi_old, coverage)
  if(coverage){  Y_new <- Eta$Etas_sample %*% t(Lambda_new) +  mvtnorm::rmvnorm(nrow(Y_old), sigma=diag(Psi_new))}
  Y_new_mean <- Eta$Etas_mean %*% t(Lambda_new)
  output <- list(Y_mean=Y_new_mean)
  if(coverage){output$Y_sample=Y_new}
  return(output)
}


predict_oos <- function(Y, Lambda, Sigma, index_impute){
  preds <- predict_y(Y[,index_impute], Lambda[-index_impute,], Lambda[index_impute,], 
                     Sigma[-index_impute], Sigma[index_impute], coverage=F)$Y_mean
  return(preds)
}


#fable_log_lik <- dmvnorm(Y_test_gene,  mu=rep(0, p), Sigma=fable_est$FABLEPostMean, log=T)
#sum(fable_log_lik)
#predict_oos(Y_test_gene, eb_init$Lambdas_eb_hat_init, fable_est$FABLEPostMean, impute_factor_index_gene)
mse_init[[g]] <- mean( (Y_test_gene[[g]][,-impute_factor_index_gene] - pred_init)^2 / sd(Y_test_gene[[g]][,-impute_factor_index_gene])^2)
print(mse_init[[g]])

eb_k_jic_log_lik <- predict_oos(Y_test_gene, eb_est_k_jic$Lambda_hat, 
                                eb_est_k_jic$sigmas_sq_mean, impute_factor_index_gene)
eb_k_jic_mse <- mean( (Y_test_gene[,-impute_factor_index_gene] - pred_init)^2 /
                        (sd(Y_test_gene[,-impute_factor_index_gene]))^2)
eb_k_jic_mse

eb_k_jic_log_lik <- predict_oos(Y_test_gene, eb_est_k_jic$Lambda_hat, eb_est_k_jic$sigmas_sq_mean, impute_factor_index_gene)
eb_k_jic_mse <- mean( (Y_test_gene[,-impute_factor_index_gene] - pred_init)^2 / (sd(Y_test_gene[,-impute_factor_index_gene]))^2)
eb_k_jic_mse

sum(eb_k_jic_log_lik)
eb_est_k_jic_over$cov_hat <- tcrossprod(eb_est_k_jic_over$Lambda_hat) + diag(eb_est_k_jic_over$sigmas_sq_mean)
eb_k_jic_over_log_lik <- dmvnorm(Y_test_gene,  mu=rep(0, p), Sigma=eb_est_k_jic_over$cov_hat, log=T)
sum(eb_k_jic_over_log_lik)

bc_est_k_jic$cov_hat <- tcrossprod(bc_est_k_jic$Lambda_hat) + diag(bc_est_k_jic$sigmas_sq_hat)
bc_k_jic_log_lik <- dmvnorm(Y_test_gene,  mu=rep(0, p), Sigma=bc_est_k_jic$cov_hat, log=T)
sum(bc_k_jic_log_lik)
bc_est_k_jic_over$cov_hat <- tcrossprod(bc_est_k_jic_over$Lambda_hat) + diag(bc_est_k_jic_over$sigmas_sq_hat)
bc_k_jic_over_log_lik <- dmvnorm(Y_test_gene,  mu=rep(0, p), Sigma=bc_est_k_jic_over$cov_hat, log=T)
sum(bc_k_jic_over_log_lik)

rot_est_k_jic$cov_hat <- tcrossprod(rot_est_k_jic$Lambda_hat) + diag(rot_est_k_jic$sigmas_sq_hat)
rot_k_jic_log_lik <- dmvnorm(Y_test_gene,  mu=rep(0, p), Sigma=rot_est_k_jic$cov_hat, log=T)
sum(rot_k_jic_log_lik)
rot_est_k_jic_over$cov_hat <- tcrossprod(rot_est_k_jic_over$Lambda_hat) + diag(rot_est_k_jic_over$sigmas_sq_hat)
rot_k_jic_over_log_lik <- dmvnorm(Y_test_gene,  mu=rep(0, p), Sigma=rot_est_k_jic_over$cov_hat, log=T)
sum(rot_k_jic_over_log_lik)

pca_est_k_jic$cov_hat <- tcrossprod(pca_est_k_jic$Lambda_hat) + diag(pca_est_k_jic$sigmas_sq_hat)
pca_k_jic_log_lik <- dmvnorm(Y_test_gene,  mu=rep(0, p), Sigma=pca_est_k_jic$cov_hat, log=T)
sum(pca_k_jic_log_lik)
pca_est_k_jic_over$cov_hat <- tcrossprod(pca_est_k_jic_over$Lambda_hat) + diag(pca_est_k_jic_over$sigmas_sq_hat)
pca_k_jic_over_log_lik <- dmvnorm(Y_test_gene,  mu=rep(0, p), Sigma=pca_est_k_jic_over$cov_hat, log=T)
sum(pca_k_jic_over_log_lik)




mse_init <- list()
mse_1 <- list()
mse_5 <- list()
mse_spectral <- list()
for(g in 1:G){
  
  print(g)
  
  pred_init <- predict_oos(Y_test_gene[[g]], eb_init$Lambdas_eb_hat_init[[g]], eb_init$Sigma_s_hat[,g], impute_factor_index_gene)
  mse_init[[g]] <- mean( (Y_test_gene[[g]][,-impute_factor_index_gene] - pred_init)^2 / sd(Y_test_gene[[g]][,-impute_factor_index_gene])^2)
  print(mse_init[[g]])
  
  pred_1 <- predict_oos(Y_test_gene[[g]], eb_one_step$Lambdas_eb_hat[[g]], eb_init$Sigma_s_hat[,g], impute_factor_index_gene)
  mse_1[[g]] <- mean( (Y_test_gene[[g]][,-impute_factor_index_gene] - pred_1)^2 / sd(Y_test_gene[[g]][,-impute_factor_index_gene])^2)
  print(mse_1[[g]])
  
  pred_5 <- predict_oos(Y_test_gene[[g]], eb_five_step$Lambdas_eb_hat[[g]], eb_init$Sigma_s_hat[,g], impute_factor_index_gene)
  mse_5[[g]] <- mean( (Y_test_gene[[g]][,-impute_factor_index_gene] - pred_5)^2 / sd(Y_test_gene[[g]][,-impute_factor_index_gene])^2)
  print(mse_5[[g]])
  
  pred_spectral <- predict_oos(Y_test_gene[[g]], eb_init$Lambdas_spectral_hat[[g]], eb_init$Sigma_s_hat[,g], impute_factor_index_gene)
  mse_spectral[[g]] <- mean( (Y_test_gene[[g]][,-impute_factor_index_gene] - pred_spectral)^2 / sd(Y_test_gene[[g]][,-impute_factor_index_gene])^2)
  print(mse_spectral[[g]])
}


eb_init <- estimate_ebss_init(Y_train_gene, k=k_s, k_max=50)
eb_one_step <- one_step_eb(
  eb_init$V_s, eb_init$D_s, eb_init$R_s, eb_init$tau_sq_s, eb_init$Sigma_s_hat, eb_init$k_0, 
  eb_init$k, update_phi = update_phi)
eb_two_step <- one_step_eb(
  eb_init$V_s, eb_init$D_s, eb_one_step$R_s, eb_one_step$tau_sq_s, eb_init$Sigma_s_hat, eb_init$k_0, 
  eb_init$k, update_phi = update_phi)

eb_three_step <- one_step_eb(
  eb_init$V_s, eb_init$D_s, eb_two_step$R_s, eb_two_step$tau_sq_s, eb_init$Sigma_s_hat, eb_init$k_0, eb_init$k
)
eb_four_step <- one_step_eb(
  eb_init$V_s, eb_init$D_s, eb_three_step$R_s, eb_three_step$tau_sq_s, eb_init$Sigma_s_hat, eb_init$k_0, eb_init$k
)
eb_five_step <- one_step_eb(
  eb_init$V_s, eb_init$D_s, eb_four_step$R_s, eb_four_step$tau_sq_s, eb_init$Sigma_s_hat, eb_init$k_0, eb_init$k
)



predict_factors <- function(Y_old, Lambda_old, Psi_old, coverage){
  k <- ncol(Lambda_old)
  Var_fact <- solve(diag(1, k, k) + t(Lambda_old) %*% diag(1/Psi_old) %*% Lambda_old)
  mean_fact <- Y_old %*% diag(1/Psi_old) %*% Lambda_old  %*% Var_fact 
  if(coverage){  eta <-mean_fact + mvtnorm::rmvnorm(nrow(Y_old), sigma=Var_fact)}
  output <- list(Etas_mean=mean_fact)
  if(coverage){output$Etas_sample=eta}
  return(output)
}

predict_y <- function(Y_old, Lambda_new, Lambda_old, Psi_new, Psi_old, coverage=F){
  Eta <- predict_factors(Y_old, Lambda_old, Psi_old, coverage)
  if(coverage){  Y_new <- Eta$Etas_sample %*% t(Lambda_new) +  mvtnorm::rmvnorm(nrow(Y_old), sigma=diag(Psi_new))}
  Y_new_mean <- Eta$Etas_mean %*% t(Lambda_new)
  output <- list(Y_mean=Y_new_mean)
  if(coverage){output$Y_sample=Y_new}
  return(output)
}


predict_oos <- function(Y, Lambda, Sigma, index_impute){
  preds <- predict_y(Y[,index_impute], Lambda[-index_impute,], Lambda[index_impute,], 
                     Sigma[-index_impute], Sigma[index_impute], coverage=F)$Y_mean
  return(preds)
}


