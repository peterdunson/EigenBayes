
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// [[Rcpp::export]]
arma::mat compute_point_estimate_Lambda_hat(
    const arma::mat& V,const arma::mat& D, const double n, const arma::vec& taus_sq,
    const arma::vec& phis_sq
) {
  const arma::uword p = V.n_rows;
  const arma::uword k = V.n_cols;
  
  if (D.n_rows != k || D.n_cols != k) {
    Rcpp::stop("D must be k x k where k = ncol(V).");
  }
  if (taus_sq.n_elem != p) {
    Rcpp::stop("taus_sq must have length p = nrow(V).");
  }
  if (phis_sq.n_elem != k) {
    Rcpp::stop("phis_sq must have length k = ncol(V).");
  }
  if (n <= 0.0) {
    Rcpp::stop("n must be positive.");
  }
  
  arma::mat Lambda_init = std::sqrt(n) * (V * D);
  
  arma::vec inv_taus = 1.0 / taus_sq;
  arma::vec inv_phis = 1.0 / phis_sq;
  arma::mat denom = n + (inv_taus * inv_phis.t());
  
  arma::mat Lambda_hat = Lambda_init / denom;
  
  return Lambda_hat;
}




// [[Rcpp::export]]
arma::cube sample_Lambda(
    const arma::mat& Lambda_hat, const arma::uvec& subsample_index, const arma::vec& taus_sq,
    const arma::vec& phis_sq, const arma::mat& sigmas_sq_samples,const double rho_sq, const double n
) {
  
  const arma::uword p_index = subsample_index.n_elem;
  const arma::uword k = Lambda_hat.n_cols;
  const arma::uword n_MC = sigmas_sq_samples.n_rows;
  
  if (sigmas_sq_samples.n_cols != p_index)
    Rcpp::stop("sigmas_sq_samples must be n_MC x p_index.");
  
  if (phis_sq.n_elem != k)
    Rcpp::stop("phis_sq must have length k.");
  
  arma::cube Z(n_MC, p_index, k);
  Z.randn();
  
  arma::cube out(n_MC, p_index, k);
  
  arma::vec inv_phis = 1.0 / phis_sq;
  arma::rowvec ones_k = arma::ones<arma::rowvec>(k);
  double sqrt_rho = std::sqrt(rho_sq);
  
  for (arma::uword t = 0; t < n_MC; ++t) {
    
    arma::vec s = sqrt_rho * arma::sqrt(sigmas_sq_samples.row(t).t()); // p_index
    
    arma::mat scale(p_index, k);
    
    for (arma::uword j = 0; j < p_index; ++j) {
      arma::uword idx = subsample_index(j)-1;
      scale.row(j) =
        s(j) / arma::sqrt(n + (1.0 / taus_sq(idx)) * inv_phis).t();
    }
    
    out.slice(t) =
      Lambda_hat.rows(subsample_index) + (Z.slice(t) % scale);
  }
  
  return out;
}


// [[Rcpp::export]]
arma::mat compute_correction_lr(
    const arma::mat& Lambda_outer, const arma::vec& sigma_sq_hat, const arma::uvec& subsample_index
  ) {
  
  const arma::uword p = Lambda_outer.n_rows;
  if (sigma_sq_hat.n_elem != p) Rcpp::stop("sigma_sq_hat must have length p.");
  
  const arma::uword p_index = subsample_index.n_elem;
  arma::mat correction(p_index, p_index, arma::fill::zeros);
  
  for (arma::uword a = 0; a < (p_index-1); ++a) {
    const arma::uword j_idx = subsample_index(a)-1;
    const double Lj = Lambda_outer(j_idx, j_idx);
    const double sj = sigma_sq_hat(j_idx);
    correction(a, a) = std::sqrt(1.0 + Lj / (2.0 * sj));
    
    for (arma::uword b = a + 1; b < p_index; ++b) {
      const arma::uword l_idx = subsample_index(b)-1;
      const double Ll = Lambda_outer(l_idx, l_idx);
      const double sj2 = sigma_sq_hat(j_idx);
      const double sl2 = sigma_sq_hat(l_idx);
      const double Ljl = Lambda_outer(j_idx, l_idx);
      
      const double num = Lj * Ll + Ljl * Ljl;
      const double den = sj2 * Ll + sl2 * Lj;
      
      const double val = std::sqrt(1.0 + num / den);
      correction(a, b) = val;
      correction(b, a) = val;
    }
  }
  
  const arma::uword j_idx = subsample_index(p_index -1)-1;
  const double Lj = Lambda_outer(j_idx, j_idx);
  const double sj = sigma_sq_hat(j_idx);
  correction((p_index -1), (p_index -1)) = std::sqrt(1.0 + Lj / (2.0 * sj));
  
  
  return correction;
}

// [[Rcpp::export]]
arma::mat compute_correction(
    const arma::mat& Lambda_outer, const arma::vec& sigma_sq_hat, const arma::uvec& subsample_index
) {
  
  const arma::uword p = Lambda_outer.n_rows;

  const arma::uword p_index = subsample_index.n_elem;
  arma::mat correction(p_index, p_index, arma::fill::zeros);
  
  for (arma::uword a = 0; a < (p_index-1); ++a) {
    const arma::uword j_idx = subsample_index(a)-1;
    const double Lj = Lambda_outer(j_idx, j_idx);
    const double sj = sigma_sq_hat(j_idx);
    correction(a, a) = std::sqrt(1.0 + Lj / (2.0 * sj));
    
    for (arma::uword b = a + 1; b < p_index; ++b) {
      const arma::uword l_idx = subsample_index(b)-1;
      const double Ll = Lambda_outer(l_idx, l_idx);
      const double sj2 = sigma_sq_hat(j_idx);
      const double sl2 = sigma_sq_hat(l_idx);
      const double Ljl = Lambda_outer(j_idx, l_idx);
      
      const double num = Lj * Ll + Ljl * Ljl;
      const double den = sj2 * Ll + sl2 * Lj;
      
      const double val = std::sqrt(1.0 + num / den);
      correction(a, b) = val;
      correction(b, a) = val;
    }
  }
  
  const arma::uword j_idx = subsample_index(p_index -1)-1;
  const double Lj = Lambda_outer(j_idx, j_idx);
  const double sj = sigma_sq_hat(j_idx);
  correction((p_index -1), (p_index -1)) = std::sqrt(1.0 + Lj / (2.0 * sj));
  
  
  return correction;
}




// [[Rcpp::export]]
arma::mat compute_sds_clt_old(
    const arma::mat& Lambda, const arma::vec& Sigma_2s
) {
  arma::mat Lambda_outer = Lambda * Lambda.t(); 
  int p = Lambda.n_rows; 
  arma::mat sds = arma::zeros<arma::mat>(p, p);
  
  for (int j = 0; j < p - 1; ++j) {
    for (int l = j + 1; l < p; ++l) {
      sds(j, l) = sqrt(Lambda_outer(j, j)*Sigma_2s(l) + Lambda_outer(l, l)*Sigma_2s(j) +
        std::pow(Lambda_outer(j, l),2) + Lambda_outer(l, l) *  Lambda_outer(j, j));
      sds(l, j) = sds(j, l); 
    }
    sds(j, j) = sqrt(2*Sigma_2s(j) + 2*std::pow(Lambda_outer(j, j),2) );
  }
  sds(p - 1, p - 1) = sqrt(2*Sigma_2s(p-1) + 
    2*std::pow(Lambda_outer(p-1, p-1),2));
  
  return sds;
}



// [[Rcpp::export]]
arma::mat compute_sds_clt_lr(const arma::mat& Lambda_outer,
                          const arma::vec& Sigma_2s) {
  
  const arma::uword p = Lambda_outer.n_rows;
  if (Sigma_2s.n_elem != p) Rcpp::stop("Sigma_2s must have length p = nrow(Lambda).");
  
  arma::mat sds(p, p, arma::fill::zeros);
  
  for (arma::uword j = 0; j + 1 < p; ++j) {
    
    for (arma::uword l = j + 1; l < p; ++l) {
      const double Ljj = Lambda_outer(j, j);
      const double Lll = Lambda_outer(l, l);
      const double Ljl = Lambda_outer(j, l);
      
      const double val =
        std::sqrt(Ljj * Sigma_2s(l) +
        Lll * Sigma_2s(j) +
        (Ljl * Ljl) +
        (Lll * Ljj));
      
      sds(j, l) = val;
      sds(l, j) = val;
    }
    
    const double Ljj = Lambda_outer(j, j);
    sds(j, j) = std::sqrt(2.0 * Sigma_2s(j) + 2.0 * (Ljj * Ljj));
  }
  
  if (p > 0) {
    const arma::uword j = p - 1;
    const double Ljj = Lambda_outer(j, j);
    sds(j, j) = std::sqrt(2.0 * Sigma_2s(j) + 2.0 * (Ljj * Ljj));
  }
  
  return sds;
}


// [[Rcpp::export]]
arma::mat compute_sds_clt(const arma::mat& Lambda_outer,
                             const arma::vec& Sigma_2s) {
  
  const arma::uword p = Lambda_outer.n_rows;
  if (Sigma_2s.n_elem != p) Rcpp::stop("Sigma_2s must have length p = nrow(Lambda).");
  
  arma::mat sds(p, p, arma::fill::zeros);
  
  for (arma::uword j = 0; j + 1 < p; ++j) {
    
    for (arma::uword l = j + 1; l < p; ++l) {
      const double Ljj = Lambda_outer(j, j);
      const double Lll = Lambda_outer(l, l);
      const double Ljl = Lambda_outer(j, l);
      
      const double val =
        std::sqrt(Ljj * Sigma_2s(l) +
        Lll * Sigma_2s(j) +
        (Ljl * Ljl) +
        (Lll * Ljj));
      
      sds(j, l) = val;
      sds(l, j) = val;
    }
    
    const double Ljj = Lambda_outer(j, j);
    //sds(j, j) = std::sqrt(2.0 * Sigma_2s(j) + 2.0 * (Ljj * Ljj));
    sds(j, j) = std::sqrt(2.0) * (Sigma_2s(j) + (Ljj * Ljj));
  }
  
  if (p > 0) {
    const arma::uword j = p - 1;
    const double Ljj = Lambda_outer(j, j);
    sds(j, j) = std::sqrt(2.0) * (Sigma_2s(j) + (Ljj * Ljj));
  }
  
  return sds;
}


// [[Rcpp::export]]
arma::mat compute_sds_clt_0(const arma::mat& Lambda_outer,
                          const double& Sigma_2) {
  
  const arma::uword p = Lambda_outer.n_rows;
  arma::mat sds(p, p, arma::fill::zeros);
  
  for (arma::uword j = 0; j + 1 < p; ++j) {
    
    for (arma::uword l = j + 1; l < p; ++l) {
      const double Ljj = Lambda_outer(j, j);
      const double Lll = Lambda_outer(l, l);
      const double Ljl = Lambda_outer(j, l);
      
      const double val =
        std::sqrt(Ljj * Sigma_2 +
        Lll * Sigma_2 +
        (Ljl * Ljl) +
        (Lll * Ljj));
      
      sds(j, l) = val;
      sds(l, j) = val;
    }
    
    const double Ljj = Lambda_outer(j, j);
    sds(j, j) = std::sqrt(2.0 * Sigma_2 + 2.0 * (Ljj * Ljj));
  }
  
  if (p > 0) {
    const arma::uword j = p - 1;
    const double Ljj = Lambda_outer(j, j);
    sds(j, j) = std::sqrt(2.0 * Sigma_2 + 2.0 * (Ljj * Ljj));
  }
  
  return sds;
}


// [[Rcpp::export]]
arma::mat compute_sds_bvm_lr(const arma::mat& Lambda_outer,
                          const arma::vec& Sigma_2s,
                          const double rho_sq = 1.0) {
  
  const arma::uword p = Lambda_outer.n_rows;
  
  if (Lambda_outer.n_cols != p)
    Rcpp::stop("Lambda_outer must be p x p.");
  if (Sigma_2s.n_elem != p)
    Rcpp::stop("Sigma_2s must have length p.");
  
  arma::mat sds(p, p, arma::fill::zeros);
  const double rho = std::sqrt(rho_sq);
  
  for (arma::uword j = 0; j + 1 < p; ++j) {
    
    const double Ljj = Lambda_outer(j, j);
    const double Sj  = Sigma_2s(j);
    
    for (arma::uword l = j + 1; l < p; ++l) {
      
      const double Lll = Lambda_outer(l, l);
      const double Sl  = Sigma_2s(l);
      
      const double val =
        rho * std::sqrt(Sl * Ljj + Sj * Lll);
      
      sds(j, l) = val;
      sds(l, j) = val;
    }
    
    sds(j, j) =
      std::sqrt(4.0 * rho_sq * Ljj * Sj + 2.0 * Sj * Sj);
  }
  
  if (p > 0) {
    const arma::uword j = p - 1;
    const double Ljj = Lambda_outer(j, j);
    const double Sj  = Sigma_2s(j);
    sds(j, j) =
      std::sqrt(4.0 * rho_sq * Ljj * Sj + 2.0 * Sj * Sj);
  }
  
  return sds;
}

// [[Rcpp::export]]
arma::mat compute_sds_bvm(const arma::mat& Lambda_outer,
                          const arma::vec& Sigma_2s,
                          const double rho_sq = 1.0) {
  
  const arma::uword p = Lambda_outer.n_rows;
  
  if (Lambda_outer.n_cols != p)
    Rcpp::stop("Lambda_outer must be p x p.");
  if (Sigma_2s.n_elem != p)
    Rcpp::stop("Sigma_2s must have length p.");
  
  arma::mat sds(p, p, arma::fill::zeros);
  const double rho = std::sqrt(rho_sq);
  
  for (arma::uword j = 0; j + 1 < p; ++j) {
    
    const double Ljj = Lambda_outer(j, j);
    const double Sj  = Sigma_2s(j);
    
    for (arma::uword l = j + 1; l < p; ++l) {
      
      const double Lll = Lambda_outer(l, l);
      const double Sl  = Sigma_2s(l);
      
      const double val =
        rho * std::sqrt(Sl * Ljj + Sj * Lll);
      
      sds(j, l) = val;
      sds(l, j) = val;
    }
    
    sds(j, j) =
      std::sqrt(4.0 * rho_sq * Ljj * Sj + 2.0 * Sj * Sj);
  }
  
  if (p > 0) {
    const arma::uword j = p - 1;
    const double Ljj = Lambda_outer(j, j);
    const double Sj  = Sigma_2s(j);
    sds(j, j) =
      std::sqrt(4.0 * rho_sq * Ljj * Sj + 2.0 * Sj * Sj);
  }
  
  return sds;
}

// [[Rcpp::export]]
arma::mat latent_factors_full_conditional_mean(
    const arma::mat& Y, const arma::mat& Lambda,const arma::vec& sigmas_sq, const arma::mat& variance
) {
  const arma::uword p = Y.n_cols;
  
  if (Lambda.n_rows != p)
    Rcpp::stop("Lambda must have nrow(Lambda)=p=ncol(Y).");
  if (sigmas_sq.n_elem != p)
    Rcpp::stop("sigmas_sq must have length p=ncol(Y).");
  
  const arma::uword k = Lambda.n_cols;
  if (variance.n_rows != k || variance.n_cols != k)
    Rcpp::stop("variance must be k x k where k=ncol(Lambda).");
  
  arma::mat Y_scaled = Y;
  Y_scaled.each_row() /= sigmas_sq.t();
  
  return (Y_scaled * Lambda) * variance;
}



// [[Rcpp::export]]
arma::mat latent_factor_mean(const arma::mat& Y,
                             const arma::cube& Lambda_samples,     // dim: n_MC x p x k
                             const arma::mat& sigmas_sq_samples) { // dim: n_MC x p
  
  const arma::uword n = Y.n_rows;
  const arma::uword p = Y.n_cols;
  
  const arma::uword n_MC = Lambda_samples.n_rows;
  const arma::uword p_ls = Lambda_samples.n_cols;
  const arma::uword k    = Lambda_samples.n_slices;
  
  if (p_ls != p) Rcpp::stop("Lambda_samples must have dimension (n_MC x p x k) with p = ncol(Y).");
  if (sigmas_sq_samples.n_rows != n_MC || sigmas_sq_samples.n_cols != p)
    Rcpp::stop("sigmas_sq_samples must be n_MC x p.");
  
  arma::mat out(n, k, arma::fill::zeros);
  arma::mat I = arma::eye<arma::mat>(k, k);
  
  for (arma::uword t = 0; t < n_MC; ++t) {
    
    arma::mat Lambda_t = Lambda_samples.row(t);
    arma::vec inv_sig = 1.0 / sigmas_sq_samples.row(t).t();
    arma::mat W = Lambda_t;
    W.each_col() %= inv_sig;
    arma::mat variance_t = arma::inv_sympd(Lambda_t.t() * W + I);
    arma::mat Y_scaled = Y;
    Y_scaled.each_row() %= inv_sig.t();
    arma::mat mean_t = (Y_scaled * Lambda_t) * variance_t;
    out += mean_t;
  }
  
  out /= static_cast<double>(n_MC);
  return out;
}


// [[Rcpp::export]]
arma::mat latent_factor_full_conditional(const arma::mat& mean_lf,
                                         const arma::mat& variance_lf) {
  
  const arma::uword n = mean_lf.n_rows;
  const arma::uword k = mean_lf.n_cols;
  
  if (variance_lf.n_rows != k || variance_lf.n_cols != k)
    Rcpp::stop("variance_lf must be k x k where k = ncol(mean_lf).");
  
  arma::mat Z = arma::randn<arma::mat>(n, k);
  arma::mat cholV = arma::chol(variance_lf); 
  
  return Z * cholV + mean_lf;
}


// [[Rcpp::export]]
Rcpp::List latent_factor_posterior_samples(
    const arma::mat& Y,
    const arma::cube& Lambda_samples,
    const arma::mat& sigmas_sq_samples
) {
  const arma::uword n = Y.n_rows;
  const arma::uword p = Y.n_cols;
  
  const arma::uword n_MC = Lambda_samples.n_rows;
  const arma::uword p_ls = Lambda_samples.n_cols;
  const arma::uword k    = Lambda_samples.n_slices;
  
  if (p_ls != p) Rcpp::stop("Lambda_samples must have dimension (n_MC x p x k) with p = ncol(Y).");
  if (sigmas_sq_samples.n_rows != n_MC || sigmas_sq_samples.n_cols != p)
    Rcpp::stop("sigmas_sq_samples must be n_MC x p.");
  
  arma::mat out_mean(n, k, arma::fill::zeros);
  
  // store samples as slices: each slice is n x k, with n_MC slices
  arma::cube out_samples(n, k, n_MC, arma::fill::zeros);
  
  arma::mat I = arma::eye<arma::mat>(k, k);
  
  for (arma::uword t = 0; t < n_MC; ++t) {
    
    // NOTE: see section (2) below about this extraction
    arma::mat Lambda_t = Lambda_samples.row(t); // likely NOT what you want dimension-wise
    
    arma::vec inv_sig = 1.0 / sigmas_sq_samples.row(t).t();
    
    arma::mat W = Lambda_t;
    W.each_col() %= inv_sig;
    
    arma::mat variance_t = arma::inv_sympd(Lambda_t.t() * W + I);
    
    arma::mat Y_scaled = Y;
    Y_scaled.each_row() %= inv_sig.t();
    
    arma::mat mean_t = (Y_scaled * Lambda_t) * variance_t;
    
    out_mean += mean_t;
    
    out_samples.slice(t) = latent_factor_full_conditional(mean_t, variance_t); // n x k
  }
  
  out_mean /= static_cast<double>(n_MC);
  
  return Rcpp::List::create(
    Rcpp::Named("posterior_mean")    = out_mean,
    Rcpp::Named("posterior_samples") = out_samples
  );
}



