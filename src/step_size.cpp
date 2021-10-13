#include "step_size.h"

// [[Rcpp::export]]
double compute_min_step_size_Beta(const List & Y_matrix_list, const List & X_list, int N) {

  R_xlen_t K = Y_matrix_list.size();

  double lipschitz = 0;

  for (R_xlen_t i = 0; i < K; i++) {

    NumericMatrix X = X_list[i];

    arma::mat X_(X.begin(), X.nrow(), X.ncol(), false);

    lipschitz += arma::accu(arma::square(X_));

  }

  NumericMatrix Y0 = Y_matrix_list[0];
  int q = Y0.ncol();

  lipschitz = std::sqrt(q) * lipschitz / N;

  return 1 / lipschitz;

}

// [[Rcpp::export]]
double compute_min_step_size_alpha(const List & Y_matrix_list) {

  NumericMatrix Y0 = Y_matrix_list[0];
  int q = Y0.ncol();

  double lipschitz = std::sqrt(q);

  return 1 / lipschitz;

}

// [[Rcpp::export]]
arma::colvec compute_min_step_size_Gamma(const List & Y_matrix_list, const List & Z_list, double rho, int N) {

  R_xlen_t K = Y_matrix_list.size();

  arma::colvec lipschitz(K);

  for (R_xlen_t i = 0; i < K; i++) {

    NumericMatrix Z = Z_list[i];

    arma::mat Z_(Z.begin(), Z.nrow(), Z.ncol(), false);

    lipschitz(i) = arma::accu(arma::square(Z_));

  }

  NumericMatrix Y0 = Y_matrix_list[0];
  int q = Y0.ncol();

  lipschitz = (std::sqrt(q) * lipschitz / N);

  return 1 / lipschitz;

}
