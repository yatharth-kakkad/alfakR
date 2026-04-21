// [[Rcpp::plugins(cpp17)]]
#include <Rcpp.h>
#include <cmath>
#include <vector>

namespace {

double log_sum_exp_cpp(const std::vector<double>& values) {
  if (values.empty()) {
    Rcpp::stop("log_sum_exp_cpp requires at least one value.");
  }
  bool has_finite_term = false;
  double max_val = R_NegInf;
  for (double value : values) {
    if (std::isnan(value) || value == R_PosInf) {
      Rcpp::stop("log_sum_exp_cpp rejects NaN and +Inf inputs.");
    }
    if (value == R_NegInf) {
      continue;
    }
    has_finite_term = true;
    if (value > max_val) {
      max_val = value;
    }
  }
  if (!has_finite_term) {
    Rcpp::stop("log_sum_exp_cpp cannot normalize an all -Inf vector.");
  }
  double accum = 0.0;
  for (double value : values) {
    if (value == R_NegInf) {
      continue;
    }
    accum += std::exp(value - max_val);
  }
  return max_val + std::log(accum);
}

double fexp_stable_cpp(double fc, double fp, double pij_value, double tt, double tol) {
  double delta = fc - fp;
  if (std::abs(delta) < tol) {
    return pij_value * fp * tt;
  }
  return pij_value * fp * std::expm1(tt * delta) / delta;
}

bool is_integer_valued_scalar(double x) {
  return std::floor(x) == x;
}

} // namespace

// [[Rcpp::export]]
Rcpp::NumericMatrix alfak_project_forward_log_cpp(Rcpp::NumericVector x0,
                                                  Rcpp::NumericVector f,
                                                  Rcpp::NumericVector timepoints) {
  const int K = x0.size();
  const int T = timepoints.size();
  if (K == 0) {
    Rcpp::stop("`x0` must contain at least one entry.");
  }
  if (f.size() != K) {
    Rcpp::stop("`x0` and `f` must have the same length.");
  }
  double x0_sum = 0.0;
  for (int i = 0; i < K; ++i) {
    if (!R_finite(x0[i]) || x0[i] < 0.0) {
      Rcpp::stop("`x0` must contain only finite non-negative values.");
    }
    if (!R_finite(f[i])) {
      Rcpp::stop("`f` must contain only finite values.");
    }
    x0_sum += x0[i];
  }
  if (!(x0_sum > 0.0) || !R_finite(x0_sum)) {
    Rcpp::stop("`x0` must sum to a positive finite value.");
  }
  for (int t = 0; t < T; ++t) {
    if (!R_finite(timepoints[t])) {
      Rcpp::stop("`timepoints` must contain only finite values.");
    }
  }
  Rcpp::NumericMatrix out(K, T);
  std::vector<double> log_x0(K);
  std::vector<double> lv(K);

  for (int i = 0; i < K; ++i) {
    log_x0[i] = std::log(x0[i] / x0_sum);
  }

  for (int t = 0; t < T; ++t) {
    for (int i = 0; i < K; ++i) {
      lv[i] = log_x0[i] + f[i] * timepoints[t];
    }
    double denom = log_sum_exp_cpp(lv);
    for (int i = 0; i < K; ++i) {
      out(i, t) = std::exp(lv[i] - denom);
    }
  }

  return out;
}

// [[Rcpp::export]]
double alfak_neg_log_lik_cpp(Rcpp::NumericVector param,
                             Rcpp::NumericMatrix counts,
                             Rcpp::NumericVector timepoints) {
  const int K = counts.nrow();
  const int T = counts.ncol();
  if (K <= 0) {
    Rcpp::stop("`counts` must have at least one row.");
  }
  if (T != timepoints.size()) {
    Rcpp::stop("`counts` must have ncol equal to length(timepoints).");
  }
  if (K == 1) {
    Rcpp::stop("`alfak_neg_log_lik_cpp()` expects at least two karyotypes; K == 1 should be handled in R.");
  }
  if (param.size() != (2 * K - 2)) {
    Rcpp::stop("`param` must have length 2*K - 2.");
  }
  std::vector<double> f_full(K, 0.0);
  std::vector<double> log_x0(K);
  std::vector<double> lv(K);

  double f_sum = 0.0;
  for (int i = 0; i < K - 1; ++i) {
    if (!R_finite(param[i])) {
      Rcpp::stop("`param` must contain only finite values.");
    }
    f_full[i] = param[i];
    f_sum += param[i];
  }
  f_full[K - 1] = -f_sum;

  for (int i = 0; i < K - 1; ++i) {
    if (!R_finite(param[K - 1 + i])) {
      Rcpp::stop("`param` must contain only finite values.");
    }
    log_x0[i] = param[K - 1 + i];
  }
  log_x0[K - 1] = 0.0;
  for (int t = 0; t < T; ++t) {
    if (!R_finite(timepoints[t])) {
      Rcpp::stop("`timepoints` must contain only finite values.");
    }
    for (int i = 0; i < K; ++i) {
      if (!R_finite(counts(i, t)) || counts(i, t) < 0.0) {
        Rcpp::stop("`counts` must contain only finite non-negative values.");
      }
    }
  }

  double nll = 0.0;
  for (int t = 0; t < T; ++t) {
    for (int i = 0; i < K; ++i) {
      lv[i] = log_x0[i] + f_full[i] * timepoints[t];
    }
    double denom = log_sum_exp_cpp(lv);
    for (int i = 0; i < K; ++i) {
      if (counts(i, t) > 0) {
        nll -= counts(i, t) * (lv[i] - denom);
      }
    }
  }

  return nll;
}

// [[Rcpp::export]]
double alfak_neighbor_objective_cpp(double fc_param,
                                    Rcpp::NumericVector parent_fitness,
                                    Rcpp::NumericVector pij_values,
                                    Rcpp::NumericVector parent_birth_times,
                                    Rcpp::NumericVector timepoints,
                                    Rcpp::NumericMatrix parent_xfit,
                                    Rcpp::NumericVector child_obs,
                                    Rcpp::NumericVector ntot,
                                    double parent_fitness_mean,
                                    double prior_mean,
                                    double prior_sd,
                                    bool do_prior,
                                    double tol) {
  const int n_parents = parent_fitness.size();
  const int n_time = timepoints.size();
  if (n_parents == 0) {
    return 1e9;
  }
  if (!R_finite(fc_param) || !R_finite(tol) || tol <= 0.0) {
    Rcpp::stop("`fc_param` must be finite and `tol` must be a positive finite value.");
  }
  if (pij_values.size() != n_parents || parent_birth_times.size() != n_parents ||
      parent_xfit.nrow() != n_parents) {
    Rcpp::stop("Parent inputs must have matching lengths/rows.");
  }
  if (parent_xfit.ncol() != n_time) {
    Rcpp::stop("`parent_xfit` must have ncol equal to length(timepoints).");
  }
  if (child_obs.size() != n_time || ntot.size() != n_time) {
    Rcpp::stop("`child_obs`, `ntot`, and `timepoints` must have matching lengths.");
  }
  if (do_prior && (!R_finite(prior_sd) || prior_sd <= 0.0 || !R_finite(prior_mean) || !R_finite(parent_fitness_mean))) {
    Rcpp::stop("When `do_prior` is TRUE, prior parameters and parent fitness mean must be finite and `prior_sd` must be positive.");
  }

  double loglik = 0.0;
  for (int t = 0; t < n_time; ++t) {
    if (!R_finite(timepoints[t])) {
      Rcpp::stop("`timepoints` must contain only finite values.");
    }
    if (!R_finite(child_obs[t]) || child_obs[t] < 0.0 || !is_integer_valued_scalar(child_obs[t])) {
      Rcpp::stop("`child_obs` must contain only finite non-negative integer-valued counts.");
    }
    if (!R_finite(ntot[t]) || ntot[t] < 0.0 || !is_integer_valued_scalar(ntot[t])) {
      Rcpp::stop("`ntot` must contain only finite non-negative integer-valued counts.");
    }
    if (child_obs[t] > ntot[t]) {
      Rcpp::stop("`child_obs` must not exceed `ntot` at any timepoint.");
    }
    double xc_est = 0.0;
    for (int p = 0; p < n_parents; ++p) {
      if (!R_finite(parent_fitness[p]) || !R_finite(pij_values[p]) || pij_values[p] < 0.0 ||
          !R_finite(parent_birth_times[p]) || !R_finite(parent_xfit(p, t))) {
        Rcpp::stop("Parent fitness, transition probabilities, birth times, and parent_xfit must be finite; pij values must be non-negative.");
      }
      double tt = std::max(0.0, timepoints[t] - parent_birth_times[p]);
      xc_est += fexp_stable_cpp(fc_param, parent_fitness[p], pij_values[p], tt, tol) * parent_xfit(p, t);
    }

    if (!R_finite(xc_est)) {
      loglik += -1e9;
      continue;
    }

    xc_est = std::max(0.0, std::min(1.0, xc_est));
    double ll = R::dbinom(child_obs[t], ntot[t], xc_est, true);
    if (!R_finite(ll)) {
      ll = -1e9;
    }
    loglik += ll;
  }

  if (do_prior && R_finite(parent_fitness_mean)) {
    double prior_ll = R::dnorm(fc_param - parent_fitness_mean, prior_mean, prior_sd, true);
    if (!R_finite(prior_ll)) {
      prior_ll = -1e9;
    }
    loglik += prior_ll;
  }

  return -loglik;
}

// [[Rcpp::export]]
Rcpp::List alfak_qr_accum_cpp(Rcpp::NumericMatrix x_trim,
                              Rcpp::NumericMatrix dx_dt) {
  const int K = x_trim.nrow();
  const int T = x_trim.ncol();
  if (K <= 0 || T < 0) {
    Rcpp::stop("`x_trim` must have positive dimensions.");
  }
  if (dx_dt.nrow() != K || dx_dt.ncol() != T) {
    Rcpp::stop("`x_trim` and `dx_dt` must have identical dimensions.");
  }
  Rcpp::NumericMatrix Q_accum(K, K);
  Rcpp::NumericVector r_accum(K);
  std::vector<double> xt(K);
  std::vector<double> xt_sq(K);
  std::vector<double> dx(K);

  for (int t = 0; t < T; ++t) {
    double sum_xt_sq = 0.0;
    double xt_dx_dot = 0.0;

    for (int i = 0; i < K; ++i) {
      if (!R_finite(x_trim(i, t)) || !R_finite(dx_dt(i, t))) {
        Rcpp::stop("`x_trim` and `dx_dt` must contain only finite values.");
      }
      xt[i] = x_trim(i, t);
      dx[i] = dx_dt(i, t);
      xt_sq[i] = xt[i] * xt[i];
      sum_xt_sq += xt_sq[i];
      xt_dx_dot += xt[i] * dx[i];
    }

    for (int i = 0; i < K; ++i) {
      r_accum[i] += xt[i] * dx[i] - xt[i] * xt_dx_dot;
      for (int j = i; j < K; ++j) {
        double value = (i == j ? xt_sq[i] : 0.0) -
          xt_sq[i] * xt[j] -
          xt[i] * xt_sq[j] +
          sum_xt_sq * xt[i] * xt[j];
        Q_accum(i, j) += value;
        if (j != i) {
          Q_accum(j, i) += value;
        }
      }
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("Q_accum") = Q_accum,
    Rcpp::Named("r_accum") = r_accum
  );
}
