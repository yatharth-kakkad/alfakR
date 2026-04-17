// [[Rcpp::plugins(cpp17)]]
#include <Rcpp.h>
#include <cmath>
#include <vector>

namespace {

double log_sum_exp_cpp(const std::vector<double>& values) {
  double max_val = values[0];
  for (double value : values) {
    if (value > max_val) {
      max_val = value;
    }
  }
  double accum = 0.0;
  for (double value : values) {
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

} // namespace

// [[Rcpp::export]]
Rcpp::NumericMatrix alfak_project_forward_log_cpp(Rcpp::NumericVector x0,
                                                  Rcpp::NumericVector f,
                                                  Rcpp::NumericVector timepoints) {
  const int K = x0.size();
  const int T = timepoints.size();
  Rcpp::NumericMatrix out(K, T);
  std::vector<double> log_x0(K);
  std::vector<double> lv(K);

  for (int i = 0; i < K; ++i) {
    log_x0[i] = std::log(x0[i]);
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
  std::vector<double> f_full(K, 0.0);
  std::vector<double> log_x0(K);
  std::vector<double> lv(K);

  if (K > 1) {
    double f_sum = 0.0;
    for (int i = 0; i < K - 1; ++i) {
      f_full[i] = param[i];
      f_sum += param[i];
    }
    f_full[K - 1] = -f_sum;
  }

  for (int i = 0; i < K; ++i) {
    log_x0[i] = param[K - 1 + i];
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

  double loglik = 0.0;
  for (int t = 0; t < n_time; ++t) {
    double xc_est = 0.0;
    for (int p = 0; p < n_parents; ++p) {
      double tt = timepoints[t] - parent_birth_times[p];
      xc_est += fexp_stable_cpp(fc_param, parent_fitness[p], pij_values[p], tt, tol) * parent_xfit(p, t);
    }

    if (!R_finite(xc_est)) {
      loglik += -1e9;
      continue;
    }

    xc_est = std::max(0.0, std::min(1.0, xc_est));
    double ll = R::dbinom(child_obs[t], std::round(ntot[t]), xc_est, true);
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
  Rcpp::NumericMatrix Q_accum(K, K);
  Rcpp::NumericVector r_accum(K);
  std::vector<double> xt(K);
  std::vector<double> xt_sq(K);
  std::vector<double> dx(K);

  for (int t = 0; t < T; ++t) {
    double sum_xt_sq = 0.0;
    double xt_dx_dot = 0.0;

    for (int i = 0; i < K; ++i) {
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
