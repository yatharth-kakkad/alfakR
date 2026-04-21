// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <sstream>
#include <cmath>
#include <numeric>
#include <limits>
using namespace Rcpp;

// [[Rcpp::interfaces(r, cpp)]]

namespace {

std::vector<int> parse_karyotype_string_cpp(const std::string& str) {
  if (str.empty()) {
    Rcpp::stop("Karyotype IDs must be non-empty character strings.");
  }
  std::vector<int> out;
  std::stringstream ss(str);
  std::string token;
  while (std::getline(ss, token, '.')) {
    if (token.empty()) {
      Rcpp::stop("Invalid karyotype ID '%s': empty chromosome-count component.", str);
    }
    std::size_t pos = 0;
    int value = 0;
    try {
      value = std::stoi(token, &pos);
    } catch (...) {
      Rcpp::stop("Invalid karyotype ID '%s': non-integer component '%s'.", str, token);
    }
    if (pos != token.size()) {
      Rcpp::stop("Invalid karyotype ID '%s': malformed component '%s'.", str, token);
    }
    if (value < 0) {
      Rcpp::stop("Invalid karyotype ID '%s': components must be non-negative integers.", str);
    }
    out.push_back(value);
  }
  if (out.empty()) {
    Rcpp::stop("Karyotype IDs must be non-empty character strings.");
  }
  return out;
}

double pij_impl(int i, int j, double beta) {
  double qij = 0.0;
  if (std::abs(i - j) > i) {
    return qij;
  }
  if (j == 0) {
    j = 2 * i;
  }
  for (int z = std::abs(i - j); z <= i; z += 2) {
    qij += R::choose(i, z) * std::pow(beta, z) * std::pow(1 - beta, i - z) *
      std::pow(0.5, z) * R::choose(z, (z + i - j) / 2);
  }
  return qij;
}

NumericMatrix validate_transition_matrix(List parms, int expected_size, const char* state_name) {
  if (!parms.containsElementNamed("A")) {
    Rcpp::stop("`parms$A` must be provided.");
  }
  SEXP a_sexp = parms["A"];
  if (!Rf_isMatrix(a_sexp) || (TYPEOF(a_sexp) != REALSXP && TYPEOF(a_sexp) != INTSXP)) {
    Rcpp::stop("`parms$A` must be a numeric matrix.");
  }
  NumericMatrix A(a_sexp);
  if (A.nrow() != A.ncol()) {
    Rcpp::stop("`parms$A` must be a square matrix.");
  }
  if (A.nrow() != expected_size) {
    Rcpp::stop("`parms$A` must have nrow(A) == ncol(A) == length(%s).", state_name);
  }
  return A;
}

} // namespace

// [[Rcpp::export]]
double pij_cpp(int i, int j, double beta) {
  if (i < 0) {
    Rcpp::stop("`i` must be a non-negative integer.");
  }
  if (j < 0) {
    Rcpp::stop("`j` must be a non-negative integer.");
  }
  if (!std::isfinite(beta) || beta < 0.0 || beta > 1.0) {
    Rcpp::stop("`beta` must be finite and in [0, 1].");
  }
  double qij = pij_impl(i, j, beta);
  if (!std::isfinite(qij) || qij < 0.0 || qij > 1.0) {
    Rcpp::stop("Internal error: computed `pij` is not finite or not in [0, 1].");
  }
  return qij;
}

 //' Prepare triplet inputs (i, j, x, dims, dimnames) for sparse A matrix.
 //' @param k_str Character vector of karyotype strings, e.g. "2.2.3".
 //' @param beta Double mis-segregation probability per chromosome.
 //' @param Nmax Optional max total mis-segregations allowed per division.
 //' If not provided, no cap is applied.
 //' @return List with elements i (rows), j (cols), x (values), dims, dimnames.
 //' @export
 // [[Rcpp::export]]
 List get_A_inputs(CharacterVector k_str, double beta, Nullable<double> Nmax_ = R_NilValue) {
   double Nmax = R_PosInf;
   if (Nmax_.isNotNull()) Nmax = as<double>(Nmax_);
   if (!std::isfinite(beta) || beta < 0.0 || beta > 1.0) {
     Rcpp::stop("`beta` must be finite and in [0, 1].");
   }
   if (!(std::isinf(Nmax) || (std::isfinite(Nmax) && Nmax >= 0.0))) {
     Rcpp::stop("`Nmax` must be Inf or a non-negative finite number.");
   }
   int n = k_str.size();
   std::vector<std::vector<int>> k_list(n);
   int num_chrom_types = -1;
   for (int i = 0; i < n; ++i) {
     std::string k_id = Rcpp::as<std::string>(k_str[i]);
     k_list[i] = parse_karyotype_string_cpp(k_id);
     if (num_chrom_types < 0) {
       num_chrom_types = static_cast<int>(k_list[i].size());
     } else if (static_cast<int>(k_list[i].size()) != num_chrom_types) {
       Rcpp::stop("All karyotype IDs must have the same number of dot-separated components.");
     }
   }
   // triplet containers (1-based indices)
   std::vector<int> ii, jj;
   std::vector<double> xx;
   std::size_t cap = static_cast<std::size_t>(n) * static_cast<std::size_t>(n);
   ii.reserve(cap);
   jj.reserve(cap);
   xx.reserve(cap);
   
   for (int i = 0; i < n; ++i) {
     const auto& ki = k_list[i];
     for (int j = 0; j < n; ++j) {
       const auto& kj = k_list[j];
       double tot = 0;
       for (size_t m = 0; m < ki.size(); ++m) tot += std::abs(ki[m] - kj[m]);
       if (tot > Nmax) continue;
       double qij = 1.0;
       for (size_t m = 0; m < ki.size(); ++m) qij *= pij_impl(ki[m], kj[m], beta);
       double val = (i == j ? (2 * qij - 1) : (2 * qij));
       if (val != 0.0) {
         ii.push_back(i + 1);
         jj.push_back(j + 1);
         xx.push_back(val);
       }
     }
   }
   // dims and dimnames
   IntegerVector dims = IntegerVector::create(n, n);
   List dimnames = List::create(k_str, k_str);
   return List::create(
     _["i"] = ii,
     _["j"] = jj,
     _["x"] = xx,
     _["dims"] = dims,
     _["dimnames"] = dimnames
   );
 }
 
// [[Rcpp::export]]
List chrmod_cpp(double time, NumericVector state, List parms) {
   NumericMatrix A = validate_transition_matrix(parms, state.size(), "state");
   int n = state.size();
   NumericVector ds(n);
   for (int j = 0; j < n; ++j) {
     double acc = 0;
     for (int i = 0; i < n; ++i) acc += state[i] * A(i, j);
     ds[j] = acc;
   }
   return List::create(ds);
 }
 
// [[Rcpp::export]]
List chrmod_rel_cpp(double time, NumericVector x, List parms) {
   NumericMatrix A = validate_transition_matrix(parms, x.size(), "x");
   int n = x.size();
   NumericVector g(n);
   for (int j = 0; j < n; ++j) {
     double acc = 0;
     for (int i = 0; i < n; ++i) acc += x[i] * A(i, j);
     g[j] = acc;
   }
   double phi = std::accumulate(g.begin(), g.end(), 0.0);
   NumericVector dx(n);
   for (int k = 0; k < n; ++k) dx[k] = g[k] - x[k] * phi;
   return List::create(dx);
 }
 
