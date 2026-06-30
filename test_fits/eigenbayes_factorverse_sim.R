# EigenBayes coverage simulation. Block loadings, k=10, n=200, p=1000.
# Runs the same 4 scenarios as the FABLE sim. H overfit to floor(sqrt(p)).
# Lorenzo's code: eigenbayes_functions.R + helper_functions_eigenbayes.cpp.

setwd("/Users/peterdunson/Desktop/EigenBayes")   # <-- adjust to repo path

source("eigenbayes_functions.R")   # sources the .cpp and defines the EB functions

out_dir   <- "test_fits"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

N         <- 200
B         <- 5
base_seed <- 42
K         <- 10
P         <- 1000
H         <- floor(sqrt(P))   # overfit upper bound = 31
alpha_ci  <- 0.05
SPHERICAL <- TRUE             # paper default; equal-Sigma sets satisfy it, unequal don't

# patch: eigenbayes_point_est uses an undefined `alpha` for the psi exponent.
# Paper uses exponent 1, so define it here in the sourced environment.
alpha <- 1

make_block_lambda <- function(P, K, val) {
  Lambda <- matrix(0, nrow = P, ncol = K)
  block  <- P / K
  for (h in seq_len(K)) {
    rows <- ((h - 1) * block + 1):(h * block)
    Lambda[rows, h] <- val
  }
  Lambda
}

sigma_unequal <- rep(c(0.1, 0.1, 0.1, 0.5, 0.5), length.out = P)
sigma_equal   <- rep(1, P)

param_sets <- list(
  list(Lambda = make_block_lambda(P, K, 1), Sigma = diag(sigma_equal)),
  list(Lambda = make_block_lambda(P, K, 1), Sigma = diag(sigma_unequal)),
  list(Lambda = make_block_lambda(P, K, 5), Sigma = diag(sigma_equal)),
  list(Lambda = make_block_lambda(P, K, 5), Sigma = diag(sigma_unequal))
)
param_sets <- lapply(param_sets, function(p) {
  p$Omega <- tcrossprod(p$Lambda) + p$Sigma
  p
})

generate_data <- function(params, N, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  L <- chol(params$Omega)
  Z <- matrix(rnorm(N * nrow(params$Lambda)), nrow = N)
  Z %*% L
}

fit_eigenbayes_fixedH <- function(dat, H, alpha_ci = 0.05, spherical = TRUE) {
  
  Y <- as.matrix(dat)
  P <- ncol(Y)
  
  t <- system.time({
    s_Y    <- svd(Y)
    eb_fit <- eigenbayes_point_est(Y, s_Y, H)
    
    # compute rho_sq on our copy so it's available for the result row
    # (eigenbayes_approx_ci computes it internally but doesn't return it)
    eb_fit <- compute_rho_sq(
      eb_fit,
      subsample_index = 1:P,
      spherical_idyosincratic_covariance = spherical
    )
    
    # credible intervals (BvM, rho-corrected) + confidence intervals (CLT)
    cis <- eigenbayes_approx_ci(
      eb_fit,
      alpha = alpha_ci,
      return_confidence_intervals = TRUE,
      return_credible_intervals   = TRUE,
      spherical_idyosincratic_covariance = spherical
    )
  })
  time_sec <- unname(t["elapsed"])
  
  LLt_hat   <- tcrossprod(eb_fit$Lambda_hat)
  Sigma_hat <- eb_fit$sigmas_sq_mean
  Omega_hat <- LLt_hat + diag(Sigma_hat)
  
  # credible-interval bounds: cis$credible_intervals is [2, p, p]
  Omega_lo <- cis$credible_intervals[1, , ]
  Omega_hi <- cis$credible_intervals[2, , ]
  
  # also keep CLT confidence intervals for comparison
  Conf_lo <- cis$confidence_intervals[1, , ]
  Conf_hi <- cis$confidence_intervals[2, , ]
  
  list(
    LLt_hat   = LLt_hat,
    Sigma_hat = Sigma_hat,
    Omega_hat = Omega_hat,
    Omega_lo  = Omega_lo,
    Omega_hi  = Omega_hi,
    Conf_lo   = Conf_lo,
    Conf_hi   = Conf_hi,
    rho_sq    = eb_fit$rho_sq,
    time_sec  = time_sec
  )
}

run_simulation_eb <- function(param_sets, N, B = 100, base_seed = 42,
                              H = 31, alpha_ci = 0.05, spherical = TRUE) {
  
  results <- vector("list", length(param_sets) * B)
  idx <- 1L
  
  for (p_idx in seq_along(param_sets)) {
    params     <- param_sets[[p_idx]]
    Omega_true <- params$Omega
    Pp         <- nrow(Omega_true)
    ut         <- upper.tri(Omega_true, diag = TRUE)
    offdiag    <- upper.tri(Omega_true, diag = FALSE)
    dg         <- diag(TRUE, Pp)
    LLt_true   <- tcrossprod(params$Lambda)
    
    for (b in seq_len(B)) {
      dat    <- generate_data(params, N, seed = base_seed + p_idx * B + b)
      result <- fit_eigenbayes_fixedH(dat, H = H, alpha_ci = alpha_ci, spherical = spherical)
      
      diff_Omega  <- result$Omega_hat - Omega_true
      diff_LLt    <- result$LLt_hat   - LLt_true
      diff_sigma2 <- result$Sigma_hat - diag(params$Sigma)
      
      LLt_hat <- result$LLt_hat
      denom   <- sqrt(sum(LLt_hat^2) * sum(LLt_true^2))
      rv_LLt  <- if (denom > 0) sum(LLt_hat * LLt_true) / denom else 0
      
      # credible-interval coverage (the rho-corrected Bayesian intervals)
      covered <- (result$Omega_lo <= Omega_true) & (Omega_true <= result$Omega_hi)
      width   <- result$Omega_hi - result$Omega_lo
      
      # CLT confidence-interval coverage, for comparison
      covered_clt <- (result$Conf_lo <= Omega_true) & (Omega_true <= result$Conf_hi)
      width_clt   <- result$Conf_hi - result$Conf_lo
      
      results[[idx]] <- data.frame(
        param_set           = p_idx, b = b,
        mse_Omega           = mean(diff_Omega^2),
        aabias_Omega        = mean(abs(diff_Omega)),
        mabias_Omega        = max(abs(diff_Omega)),
        mse_LLt             = mean(diff_LLt^2),
        aabias_LLt          = mean(abs(diff_LLt)),
        mabias_LLt          = max(abs(diff_LLt)),
        mse_sigma2          = mean(diff_sigma2^2),
        aabias_sigma2       = mean(abs(diff_sigma2)),
        mabias_sigma2       = max(abs(diff_sigma2)),
        rv_LLt              = rv_LLt,
        cover_Omega         = mean(covered[ut]),
        cover_Omega_diag    = mean(covered[dg]),
        cover_Omega_offdiag = if (any(offdiag)) mean(covered[offdiag]) else NA,
        ci_width_Omega      = mean(width[ut]),
        cover_clt           = mean(covered_clt[ut]),
        ci_width_clt        = mean(width_clt[ut]),
        rho_sq              = result$rho_sq,
        time_sec            = result$time_sec
      )
      idx <- idx + 1L
    }
  }
  
  do.call(rbind, results)
}

results_eb <- run_simulation_eb(
  param_sets = param_sets,
  N          = N,
  B          = B,
  base_seed  = base_seed,
  H          = H,
  alpha_ci   = alpha_ci,
  spherical  = SPHERICAL
)

results_eb$method <- "eigenbayes"

fname <- file.path(out_dir, sprintf("sim_eigenbayes_%d_%d.rds", K, P))
saveRDS(results_eb, fname)

set_label <- c("lambda=1, Sigma=I", "lambda=1, Sigma uneq",
               "lambda=5, Sigma=I", "lambda=5, Sigma uneq")

for (ps in sort(unique(results_eb$param_set))) {
  r <- results_eb[results_eb$param_set == ps, ]
  print(data.frame(
    set        = set_label[ps],
    cover_cred = round(mean(r$cover_Omega), 3),
    cover_clt  = round(mean(r$cover_clt), 3),
    ci_width   = round(mean(r$ci_width_Omega), 3),
    rho2       = round(mean(r$rho_sq), 2),
    mse_Omega  = round(mean(r$mse_Omega), 4)
  ), row.names = FALSE)
}
