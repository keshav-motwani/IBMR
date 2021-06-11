#include "probabilities.h"
#include "gradient.h"

// [[Rcpp::export]]
arma::mat compute_gradient_Beta(const List & Y_matrix_list, const List & X_list, const List & Z_list, const arma::colvec & alpha, const arma::mat & Beta, const std::vector<arma::mat> & Gamma_list, int N) {

  R_xlen_t K = Y_matrix_list.size();

  arma::mat gradient(Beta.n_rows, Beta.n_cols);
  gradient.zeros();

  for (R_xlen_t i = 0; i < K; i++) {

    NumericMatrix Y = Y_matrix_list[i];
    NumericMatrix X = X_list[i];
    NumericMatrix Z = Z_list[i];
    const arma::mat & Gamma = Gamma_list[i];

    arma::mat Y_(Y.begin(), Y.nrow(), Y.ncol(), false);
    arma::mat X_(X.begin(), X.nrow(), X.ncol(), false);
    arma::mat Z_(Z.begin(), Z.nrow(), Z.ncol(), false);

    arma::mat P = compute_probabilities(X_, Z_, alpha, Beta, Gamma);
    arma::mat C = compute_conditional_probabilities(Y_, P);

    gradient += X_.t() * (P - C);

  }

  return gradient / N;

}

// [[Rcpp::export]]
arma::colvec compute_gradient_alpha(const List & Y_matrix_list, const List & X_list, const List & Z_list, const arma::colvec & alpha, const arma::mat & Beta, const std::vector<arma::mat> & Gamma_list, int N) {

  R_xlen_t K = Y_matrix_list.size();

  arma::colvec gradient(Beta.n_cols);
  gradient.zeros();

  for (R_xlen_t i = 0; i < K; i++) {

    NumericMatrix Y = Y_matrix_list[i];
    NumericMatrix X = X_list[i];
    NumericMatrix Z = Z_list[i];
    const arma::mat & Gamma = Gamma_list[i];

    arma::mat Y_(Y.begin(), Y.nrow(), Y.ncol(), false);
    arma::mat X_(X.begin(), X.nrow(), X.ncol(), false);
    arma::mat Z_(Z.begin(), Z.nrow(), Z.ncol(), false);

    arma::mat P = compute_probabilities(X_, Z_, alpha, Beta, Gamma);
    arma::mat C = compute_conditional_probabilities(Y_, P);

    gradient += arma::sum(P - C, 0).t();

  }

  return gradient / N;

}

// [[Rcpp::export]]
arma::mat compute_gradient_Gamma(const arma::mat & Y, const arma::mat & X, const arma::mat & Z, const arma::colvec & alpha, const arma::mat & Beta, const arma::mat & Gamma, double rho, int N) {

  arma::mat P = compute_probabilities(X, Z, alpha, Beta, Gamma);
  arma::mat C = compute_conditional_probabilities(Y, P);

  return Z.t() * (P - C) / N + rho * Gamma;

}
