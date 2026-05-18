

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
bc_k_jic <- data.frame()
bc_k_jic_over <- data.frame()
bc_k_over <- data.frame()
bc_k_sqrt_p <- data.frame()

eb_true_k <- data.frame()
eb_k_jic <- data.frame()
eb_k_jic_over <- data.frame()
eb_k_over <- data.frame()
eb_k_sqrt_p <- data.frame()

rotate_true_k <- data.frame()
rotate_k_jic <- data.frame()
rotate_k_jic_over <- data.frame()
rotate_k_over <- data.frame()
rotate_k_sqrt_p <- data.frame()

pca_true_k <- data.frame()
pca_k_jic <- data.frame()
pca_k_jic_over <- data.frame()
pca_k_over <- data.frame()
pca_k_sqrt_p <- data.frame()

fable_k_hat <- data.frame()
test_barigozzi_cho <- T; test_eigenbayes <- T; test_fable <- T; test_pca <- T; test_rotate <- T

adapt_to_outcome <- F; sparse <- T

if(adapt_to_outcome){
  scenario <- scenario + 4
}
if(sparse){
  scenario <- scenario + 8
}

k_sqrt_p <- floor(sqrt(p))


n_sim <- 50
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
  k_over <- floor(sqrt(n))
  
  if(test_fable){
    ptm <- proc.time() 
    fable_est <- FABLEPosteriorMean_2(data$Y)
    fable_time <- proc.time() - ptm
    fable_est$estRank
    fable_k_hat <- rbind(fable_k_hat, c(compute_metrics_fable(fable_est, compute_coverage_=T, idx_cvg=idx_cvg, idx_fr=idx_fr),
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
    # k_hat jic
    ptm <- proc.time() 
    bc_est <- barigozzi_cho_est(data$Y, k_jic)
    bc_time <- proc.time() - ptm
    bc_k_jic <- rbind(bc_k_jic, c(compute_metrics(bc_est, idx_cvg=idx_cvg, idx_fr=idx_fr), bc_time[3], k_jic))
    # k_hat jic + 5
    ptm <- proc.time() 
    bc_est <- barigozzi_cho_est(data$Y, k_jic + 10)
    bc_time <- proc.time() - ptm
    bc_k_jic_over <- rbind(bc_k_jic_over, c(compute_metrics(bc_est, idx_cvg=idx_cvg, idx_fr=idx_fr), bc_time[3], k_jic +10))
    # k + 5
    ptm <- proc.time() 
    bc_est <- barigozzi_cho_est(data$Y, k_over)
    bc_time <- proc.time() - ptm
    bc_k_over <- rbind(bc_k_over, c(compute_metrics(bc_est, idx_cvg=idx_cvg, idx_fr=idx_fr), bc_time[3], k+5))
    rm(bc_est)
    # k sqrt p
    ptm <- proc.time() 
    bc_est <- barigozzi_cho_est(data$Y, k_sqrt_p)
    bc_time <- proc.time() - ptm
    bc_k_sqrt_p <- rbind(bc_k_sqrt_p, c(compute_metrics(bc_est, idx_cvg=idx_cvg, idx_fr=idx_fr), bc_time[3], k_sqrt_p))
    
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
    eb_true_k <- rbind(eb_true_k, c(compute_metrics(eb_est, compute_coverage_=T, idx_cvg=idx_cvg, idx_fr=idx_fr), 
                                    eb_time[3], k))
    eb_true_k
    # k_hat jic
    #ptm <- proc.time() 
    eb_est <- eigenbayes_point_est(data$Y, s_Y, k_jic)
    #eb_time <- proc.time() - ptm
    eb_k_jic <- rbind(eb_k_jic, c(compute_metrics(eb_est, compute_coverage_=T, idx_cvg=idx_cvg, idx_fr=idx_fr), 
                                  eb_time[3], k_jic))
    # k_hat jic + 5
    #ptm <- proc.time() 
    eb_est <- eigenbayes_point_est(data$Y, s_Y, k_jic+10)
    #eb_time <- proc.time() - ptm
    eb_k_jic_over <- rbind(eb_k_jic_over, c(compute_metrics(eb_est, compute_coverage_=T, idx_cvg=idx_cvg, idx_fr=idx_fr),
                                            eb_time[3], k_jic +10))
    # k + 5
    #ptm <- proc.time() 
    eb_est <- eigenbayes_point_est(data$Y, s_Y, k_over)
    #eb_time <- proc.time() - ptm
    eb_k_over <- rbind(eb_k_over, c(compute_metrics(eb_est, compute_coverage_=T, idx_cvg=idx_cvg, idx_fr=idx_fr),
                                    eb_time[3], k+5))
    rm(eb_est)
    # k sqrt p
    #ptm <- proc.time() 
    eb_est <- eigenbayes_point_est(data$Y, s_Y, k_sqrt_p)
    #eb_time <- proc.time() - ptm
    eb_k_sqrt_p <- rbind(eb_k_sqrt_p, c(compute_metrics(eb_est, compute_coverage_=T, idx_cvg=idx_cvg, idx_fr=idx_fr),
                                        eb_time[3], k_sqrt_p))
    
  }
  
  if(test_rotate){
    # true k
    ptm <- proc.time() 
    rot_est <- rotate_est(data$Y, k)
    rotate_time <- proc.time() - ptm
    rotate_true_k <- rbind(rotate_true_k, c(compute_metrics(rot_est, idx_cvg=idx_cvg, idx_fr=idx_fr), rotate_time[3], k))
    # k_hat jic
    ptm <- proc.time() 
    rot_est <- rotate_est(data$Y, k_jic)
    rot_time <- proc.time() - ptm
    rotate_k_jic <- rbind(rotate_k_jic, c(compute_metrics(rot_est, idx_cvg=idx_cvg, idx_fr=idx_fr), rot_time[3], k_jic))
    # k_hat jic + 10
    ptm <- proc.time() 
    rot_est <- rotate_est(data$Y, k_jic + 10)
    rot_time <- proc.time() - ptm
    rotate_k_jic_over <- rbind(rotate_k_jic_over, c(compute_metrics(rot_est, idx_cvg=idx_cvg, idx_fr=idx_fr), rot_time[3], k_jic + 10))
    # k + 5
    ptm <- proc.time() 
    rot_est <- rotate_est(data$Y, k_over)
    rotate_time <- proc.time() - ptm
    rotate_k_over <- rbind(rotate_k_over, c(compute_metrics(rot_est,idx_cvg=idx_cvg, idx_fr=idx_fr), rotate_time[3], k+5))
    rm(rot_est)
    # k sqrt p
    ptm <- proc.time() 
    rot_est <- rotate_est(data$Y, k_sqrt_p)
    rotate_time <- proc.time() - ptm
    rotate_k_sqrt_p <- rbind(rotate_k_sqrt_p, c(compute_metrics(rot_est,idx_cvg=idx_cvg, idx_fr=idx_fr), rotate_time[3], k_sqrt_p))
  }
  
  if(test_pca){
    # true k
    ptm <- proc.time() 
    pca_est <- spectral_est(data$Y, s_Y, k)
    pca_time <- proc.time() - ptm
    pca_true_k <- rbind(pca_true_k, c(compute_metrics(pca_est, idx_cvg=idx_cvg, idx_fr=idx_fr), pca_time[3], k))
    # k_hat jic
    ptm <- proc.time() 
    pca_est <- spectral_est(data$Y, s_Y, k_jic)
    pca_time <- proc.time() - ptm
    pca_k_jic <- rbind(pca_k_jic, c(compute_metrics(pca_est, idx_cvg=idx_cvg, idx_fr=idx_fr), pca_time[3], k_jic))
    # k_hat jic + 10
    ptm <- proc.time() 
    pca_est <- spectral_est(data$Y, s_Y, k_jic + 10)
    pca_time <- proc.time() - ptm
    pca_k_jic_over <- rbind(pca_k_jic_over, c(compute_metrics(pca_est, idx_cvg=idx_cvg, idx_fr=idx_fr), pca_time[3], k_jic + 10))
    # k + 5
    ptm <- proc.time() 
    pca_est <- spectral_est(data$Y, s_Y, k_over)
    pca_time <- proc.time() - ptm
    pca_k_over <- rbind(pca_k_over, c(compute_metrics(pca_est, idx_cvg=idx_cvg, idx_fr=idx_fr), pca_time[3], k_over))
    # k sqrt p
    ptm <- proc.time() 
    pca_est <- spectral_est(data$Y, s_Y, k_sqrt_p)
    pca_time <- proc.time() - ptm
    pca_k_sqrt_p <- rbind(pca_k_sqrt_p, c(compute_metrics(pca_est, idx_cvg=idx_cvg, idx_fr=idx_fr), pca_time[3], k_sqrt_p))
  }
  
}




names <- c('L_fr', 'L_fr_sub',  'C_fr', 'C_fr_sub', 'X_fr', 'X_fr_sub', 'ci_cov', 'ci_len', 'cr_cov', 'cr_len', 'time', 'k')

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

names(bc_k_sqrt_p)<- names
names(eb_k_sqrt_p)<- names
names(rotate_k_sqrt_p)<- names
names(pca_k_sqrt_p)<- names



names(fable_k_hat) <- names

dir.create(file.path("simulations", scenario), recursive = TRUE, showWarnings = FALSE)

write.csv(bc_true_k, file.path("simulations", scenario, "bc_true_k.csv"))
write.csv(bc_k_jic, file.path("simulations", scenario, "bc_k_jic.csv"))
write.csv(bc_k_jic_over, file.path("simulations", scenario, "bc_k_jic_over.csv"))
write.csv(bc_k_over, file.path("simulations", scenario, "bc_k_over.csv"))

write.csv(eb_true_k, file.path("simulations", scenario, "eb_true_k.csv"))
write.csv(eb_k_jic, file.path("simulations", scenario, "eb_k_jic.csv"))
write.csv(eb_k_jic_over, file.path("simulations", scenario, "eb_k_jic_over.csv"))
write.csv(eb_k_over, file.path("simulations", scenario, "eb_k_over.csv"))

write.csv(rotate_true_k, file.path("simulations", scenario, "rotate_true_k.csv"))
write.csv(rotate_k_jic, file.path("simulations", scenario, "rotate_k_jic.csv"))
write.csv(rotate_k_jic_over, file.path("simulations", scenario, "rotate_k_jic_over.csv"))
write.csv(rotate_k_over, file.path("simulations", scenario, "rotate_k_over.csv"))

write.csv(pca_true_k, file.path("simulations", scenario, "pca_true_k.csv"))
write.csv(pca_k_jic, file.path("simulations", scenario, "pca_k_jic.csv"))
write.csv(pca_k_jic_over, file.path("simulations", scenario, "pca_k_jic_over.csv"))
write.csv(pca_k_over, file.path("simulations", scenario, "pca_k_over.csv"))

write.csv(bc_k_sqrt_p, file.path("simulations", scenario, "bc_k_sqrt_p.csv"))
write.csv(eb_k_sqrt_p, file.path("simulations", scenario, "eb_k_sqrt_p.csv"))
write.csv(rotate_k_sqrt_p, file.path("simulations", scenario, "rotate_k_sqrt_p.csv"))
write.csv(pca_k_sqrt_p, file.path("simulations", scenario, "pca_k_sqrt_p.csv"))

write.csv(fable_k_hat, file.path("simulations", scenario, "fable_k_hat.csv"))







library(ggplot2)

#scenario <- 1
sc_dir <- file.path("simulations", scenario)

read_one_res <- function(path){
  tmp <- read.csv(path, check.names = FALSE)
  vals <- tmp$C_fr
  vals <- vals[!is.na(vals)]
  
  file <- tools::file_path_sans_ext(basename(path))
  
  method <- if(grepl("^bc_", file)) {
    "BC"
  } else if(grepl("^eb_", file)) {
    "EB"
  } else if(grepl("^rotate_", file)) {
    "ROTATE"
  } else if(grepl("^pca_", file)) {
    "PCA"
  } else if(grepl("^fable_", file)) {
    "FABLE"
  }
  
  variant <- sub("^[^_]+_", "", file)
  variant_lab <- switch(
    variant,
    "true_k"     = "k",
    "k_jic_over" = "JIC+10",
    "k_sqrt_p"   = "sqrt(p)",
    "k_hat"      = "JIC",
    variant
  )
  
  label <- if(method == "FABLE") "FABLE" else paste0(method, " (", variant_lab, ")")
  
  data.frame(
    value = vals,
    Method = method,
    label = label
  )
}

label_keep <- c(
  "EB (k)", "EB (JIC+10)", "EB (sqrt(p))",
  "BC (k)", "BC (JIC+10)", "BC (sqrt(p))",
  "PCA (k)", "PCA (JIC+10)", "PCA (sqrt(p))",
  "ROTATE (k)", "ROTATE (JIC+10)", "ROTATE (sqrt(p))",
  "FABLE"
)

csv_files <- list.files(sc_dir, pattern = "\\.csv$", full.names = TRUE)
plot_dat <- do.call(rbind, lapply(csv_files, read_one_res))

plot_dat <- subset(plot_dat, label %in% label_keep)

plot_dat$Method <- factor(plot_dat$Method, levels = c("EB", "BC", "PCA", "ROTATE", "FABLE"))
plot_dat$label  <- factor(plot_dat$label, levels = label_keep[label_keep %in% unique(plot_dat$label)])

p1 <- ggplot(plot_dat, aes(x = label, y = value, fill = Method)) +
  geom_boxplot(outlier.size = 0.6) +
  theme_bw() +
  labs(x = "Method", y = "Loss") +
  scale_x_discrete(labels = c(
    "BC (k)" = "BC (k)",
    "BC (JIC+10)" = "BC (JIC+10)",
    "BC (sqrt(p))" = expression(BC~"("~sqrt(p)~")"),
    "EB (k)" = "EB (k)",
    "EB (JIC+10)" = "EB (JIC+10)",
    "EB (sqrt(p))" = expression(EB~"("~sqrt(p)~")"),
    "ROTATE (k)" = "ROTATE (k)",
    "ROTATE (JIC+10)" = "ROTATE (JIC+10)",
    "ROTATE (sqrt(p))" = expression(ROTATE~"("~sqrt(p)~")"),
    "PCA (k)" = "PCA (k)",
    "PCA (JIC+10)" = "PCA (JIC+10)",
    "PCA (sqrt(p))" = expression(PCA~"("~sqrt(p)~")"),
    "FABLE" = "FABLE"
  )) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8)
  )

p1
 

ggsave(
  file.path('simulations/fig/', paste0("scenario_", scenario, ".png")),
  plot = p1,
  width = 16,
  height = 8,
  dpi = 300
)

p1

read_one_res <- function(path){
  tmp <- read.csv(path, check.names = FALSE)
  vals <- tmp[[ncol(tmp)]]
  vals <- vals[!is.na(vals)]
  
  file <- tools::file_path_sans_ext(basename(path))
  
  method <- if(grepl("^bc_", file)) {
    "BC"
  } else if(grepl("^eb_", file)) {
    "EB"
  } else if(grepl("^rotate_", file)) {
    "ROTATE"
  } else if(grepl("^pca_", file)) {
    "PCA"
  } else if(grepl("^fable_", file)) {
    "FABLE"
  }
  
  variant <- sub("^[^_]+_", "", file)
  variant_lab <- switch(
    variant,
    "true_k"     = "k",
    "k_jic_over" = "JIC+10",
    "k_sqrt_p"   = "sqrt(p)",
    "k_hat"      = "k hat",
    variant
  )
  
  label <- if(method == "FABLE") "FABLE" else paste0(method, " (", variant_lab, ")")
  
  data.frame(
    value = vals,
    Method = method,
    label = label
  )
}

label_keep <- c(
  "BC (k)", "BC (JIC+10)", "BC (sqrt(p))",
  "EB (k)", "EB (JIC+10)", "EB (sqrt(p))",
  "ROTATE (k)", "ROTATE (JIC+10)", "ROTATE (sqrt(p))",
  "PCA (k)", "PCA (JIC+10)", "PCA (sqrt(p))",
  "FABLE"
)

scenario_dirs <- list.dirs("simulations", recursive = FALSE, full.names = TRUE)

for(sc_dir in scenario_dirs){
  
  csv_files <- list.files(sc_dir, pattern = "\\.csv$", full.names = TRUE)
  plot_dat <- do.call(rbind, lapply(csv_files, read_one_res))
  
  plot_dat <- subset(plot_dat, label %in% label_keep)
  
  plot_dat$Method <- factor(plot_dat$Method, levels = c("BC", "EB", "PCA", "ROTATE", "FABLE"))
  plot_dat$label  <- factor(plot_dat$label, levels = label_keep[label_keep %in% unique(plot_dat$label)])
  
  p <- ggplot(plot_dat, aes(x = label, y = value, fill = Method)) +
    geom_boxplot(outlier.size = 0.6) +
    theme_bw() +
    labs(x = "Method", y = "") +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 8)
    )
  
  ggsave(
    file.path(sc_dir, paste0(basename(sc_dir), "_boxplot.png")),
    plot = p,
    width = 10,
    height = 4.5,
    dpi = 300
  )
}




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
