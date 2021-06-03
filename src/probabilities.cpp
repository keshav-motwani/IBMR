#include "probabilities.h"

//' @export
// [[Rcpp::export]]
arma::mat compute_probabilities(const arma::mat & X, const arma::mat & Z, const arma::colvec & alpha, const arma::mat & Beta, const arma::mat & Gamma) {

  arma::colvec o = arma::ones<arma::colvec>(X.n_rows);

  arma::mat P = arma::exp(o * alpha.t() + X * Beta + Z * Gamma);

  return P.each_col() % (1 / arma::sum(P, 1));

}

//' @export
// [[Rcpp::export]]
arma::mat compute_conditional_probabilities(const arma::mat & Y_matrix, const arma::mat & P) {

  arma::mat C = P % Y_matrix;

  return C.each_col() % (1 / arma::sum(C, 1));

}
