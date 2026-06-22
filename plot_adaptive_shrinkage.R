k <- 10
p <- 5000
n <- 1000
data <- gen_factor_data(
  n=n, p=p, psi_min=1, psi_max=5, r=k, 
  loading_scales=seq(0.15, 0.25, length.out = k),
  adapt_to_outcome = T, sparse = F
)


names(data)
X <- data$F_mat %*% t(data$Lambda)
sparse_idx <- data$sparse_idx
Lambda_outer_0 <- tcrossprod(data$Lambda)
Theta_0 <- Lambda_outer_0 + diag(data$Sigma)
boxplot(diag(Lambda_outer_0) / data$Sigma)

s_Lambda <- svd(data$Lambda)
plot(s_Lambda$d)

hypers_save_1 <- matrix(NA, 50, k+5)
hypers_save_2 <- matrix(NA, 50, p)

for(t in 1:50){
  Y <- matrix(rnorm(n*k), ncol=k) %*% t(data$Lambda) + matrix(rnorm(n*p), ncol=p) %*% diag(sqrt(data$Sigma))
  s_Y <- svd(Y)
  eb_est <- eigenbayes_point_est(data$Y, s_Y, k+5)
  #plot(sum(eb_est$taus_sq * eb_est$sigmas_sq_hat) * eb_est$phis_sq)
  #points(c(s_Lambda$d^2, rep(0, 5)), col='red')
  hypers_save_1[t,] <- sum(eb_est$taus_sq * eb_est$sigmas_sq_hat) * eb_est$phis_sq
  hypers_save_2[t,] <- (eb_est$taus_sq * eb_est$sigmas_sq_hat) * sum(eb_est$phis_sq)
  
  rm(eb_est, s_Y)
}


cols_norm_sq <- rowSums(data$Lambda^2)

plot(cols_norm_sq , hypers_save_2[1,]); abline(0,1, col='red')


library(ggplot2)
library(tidyr)
library(dplyr)

# reshape data
df_long <- as.data.frame(hypers_save_1) %>%
  setNames(as.character(1:ncol(hypers_save_1))) %>%
  pivot_longer(cols = everything(), names_to = "component", values_to = "value") %>%
  mutate(component = factor(component, levels = as.character(1:ncol(hypers_save_1))))

# true values
vals <- c(s_Lambda$d^2, rep(0, 5))

df_points <- data.frame(
  component = factor(1:length(vals), levels = 1:ncol(hypers_save_1)),
  value = vals,
  type = "True parameter value"
)

df_emp_est <- data.frame(
  component = factor(1:length(vals), levels = 1:ncol(hypers_save_1)),
  value = s_Y$d[1:(k+5)]^2 / n,
  type = "Empirical Estimates"
)

df_long$type <- "Expected value a priori"

ggplot() +
  geom_boxplot(
    data = df_long,
    aes(x = component, y = value, fill = type)
  ) +
  geom_point(
    data = df_points,
    aes(x = component, y = value, color = type),
    size = 2
  ) +
  scale_fill_manual(values = c("Expected value a priori" = "grey75")) +
  scale_color_manual(values = c("True parameter value" = "red")) +
  labs(
    title = "",
    x = "Latent component",
    y = "Magnitude",
    fill = "",
    color = ""
  ) +
  theme_minimal(base_size = 18) +   # <-- global scaling
  theme(
    plot.title = element_text(hjust = 0.5, size = 24),
    axis.title = element_text(size = 22),
    axis.text = element_text(size = 18),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 18),
    legend.position = "top"
  )

ggsave(
  filename = "fig/adaptive_shrinkage_latent.png",
  plot = last_plot(),
  width = 8.27,
  height = 5.85,
  units = "in",
  dpi = 600
)


df_plot <- data.frame(
  true_value = cols_norm_sq,
  estimate = hypers_save_2[1,]
)

ggplot(df_plot, aes(x = true_value, y = estimate)) +
  geom_point(size = 2) +
  geom_abline(
    intercept = 0,
    slope = 1,
    color = "red",
    linewidth = 1
  ) +
  labs(
    title = "",
    x = "True parameter value",
    y = "Expected value a priori"
  ) +
  theme_minimal(base_size = 18) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 24),
    axis.title = element_text(size = 22),
    axis.text = element_text(size = 18),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 18),
    legend.position = "top"
  )

ggsave(
  filename = "fig/adaptive_shrinkage_outcome.png",
  plot = last_plot(),
  width = 8.27,
  height = 5.85,
  units = "in",
  dpi = 600
)



mean(diag(Lambda_outer_0) / data$Sigma )
sum(X^2) / sum((data$Y - X)^2)
