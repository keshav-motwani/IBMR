#include "probabilities.h"
#include "objective.h"

// [[Rcpp::export]]
double compute_negative_log_likelihood(const List & Y_matrix_list, const List & X_list, const List & Z_list, const arma::colvec & alpha, const arma::mat & Beta, const arma::field<arma::mat> & Gamma_list, int N) {

  R_xlen_t K = Y_matrix_list.size();

  double ll = 0;

  for (R_xlen_t i = 0; i < K; i++) {

    NumericMatrix Y = Y_matrix_list[i];
    NumericMatrix X = X_list[i];
    NumericMatrix Z = Z_list[i];
    const arma::mat & Gamma = Gamma_list(i);

    arma::mat Y_(Y.begin(), Y.nrow(), Y.ncol(), false);
    arma::mat X_(X.begin(), X.nrow(), X.ncol(), false);
    arma::mat Z_(Z.begin(), Z.nrow(), Z.ncol(), false);

    arma::mat P = arma::exp(compute_linear_predictor(X_, Z_, alpha, Beta, Gamma));

    ll += arma::accu(arma::log(arma::sum(P % Y_, 1)) - arma::log(arma::sum(P, 1)));

  }

  return -1 * ll / N;

}

//' @export
// [[Rcpp::export]]
double compute_negative_log_likelihood_no_Gamma(const List & Y_matrix_list, const List & X_list, const arma::colvec & alpha, const arma::mat & Beta, int N) {

  R_xlen_t K = Y_matrix_list.size();

  double ll = 0;

  for (R_xlen_t i = 0; i < K; i++) {

    NumericMatrix Y = Y_matrix_list[i];
    NumericMatrix X = X_list[i];

    arma::mat Y_(Y.begin(), Y.nrow(), Y.ncol(), false);
    arma::mat X_(X.begin(), X.nrow(), X.ncol(), false);

    arma::mat P = arma::exp(compute_linear_predictor_no_Gamma(X_, alpha, Beta));

    ll += arma::accu(arma::log(arma::sum(P % Y_, 1)) - arma::log(arma::sum(P, 1)));

  }

  return -1 * ll / N;

}

//' @export
// [[Rcpp::export]]
double compute_negative_log_likelihood_from_probabilities(const List & Y_matrix_list, const List & P_list, int N) {

  R_xlen_t K = Y_matrix_list.size();

  double ll = 0;

  for (R_xlen_t i = 0; i < K; i++) {

    NumericMatrix Y = Y_matrix_list[i];
    NumericMatrix P = P_list[i];

    arma::mat Y_(Y.begin(), Y.nrow(), Y.ncol(), false);
    arma::mat P_(P.begin(), P.nrow(), P.ncol(), false);

    ll += arma::accu(arma::log(arma::sum(P_ % Y_, 1)));

  }

  return -1 * ll / N;

}

// [[Rcpp::export]]
double compute_negative_log_likelihood_1(const arma::mat & Y, const arma::mat & X, const arma::mat & Z, const arma::colvec & alpha, const arma::mat & Beta, const arma::mat & Gamma, int N) {

  arma::mat P = arma::exp(compute_linear_predictor(X, Z, alpha, Beta, Gamma));

  return -1 * arma::accu(arma::log(arma::sum(P % Y, 1)) - arma::log(arma::sum(P, 1))) / N;

}

// [[Rcpp::export]]
double group_lasso_penalty(const arma::mat & Beta, double lambda) {

  return lambda * arma::accu(arma::sqrt(arma::sum(arma::square(Beta), 1)));

}

// [[Rcpp::export]]
double l2_penalty(const arma::field<arma::mat> & Gamma_list, double rho) {

  R_xlen_t K = Gamma_list.size();

  double penalty = 0;

  for (R_xlen_t i = 0; i < K; i++) {

    penalty += arma::accu(arma::square(Gamma_list(i)));

  }

  return (rho / 2) * penalty;

}

// [[Rcpp::export]]
double compute_objective_function(const List & Y_matrix_list, const List & X_list, const List & Z_list, const arma::colvec & alpha, const arma::mat & Beta, const arma::field<arma::mat> & Gamma_list, double lambda, double rho, int N) {

  double nll = compute_negative_log_likelihood(Y_matrix_list, X_list, Z_list, alpha, Beta, Gamma_list, N);
  double gl = group_lasso_penalty(Beta, lambda);
  double l2 = l2_penalty(Gamma_list, rho);

  return nll + gl + l2;

}
