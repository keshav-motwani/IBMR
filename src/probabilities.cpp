#include "probabilities.h"

// [[Rcpp::export]]
arma::mat compute_linear_predictor(const arma::mat & X, const arma::mat & Z, const arma::colvec & alpha, const arma::mat & Beta, const arma::mat & Gamma) {

  arma::mat lin_pred = X * Beta + Z * Gamma;
  lin_pred = lin_pred.each_row() + alpha.t();

  return lin_pred;

}

// [[Rcpp::export]]
arma::mat compute_linear_predictor_no_Gamma(const arma::mat & X, const arma::colvec & alpha, const arma::mat & Beta) {

  arma::mat lin_pred = X * Beta;
  lin_pred = lin_pred.each_row() + alpha.t();

  return lin_pred;

}

// [[Rcpp::export]]
arma::mat compute_probabilities(const arma::mat & X, const arma::mat & Z, const arma::colvec & alpha, const arma::mat & Beta, const arma::mat & Gamma) {

  arma::mat P = arma::exp(compute_linear_predictor(X, Z, alpha, Beta, Gamma));

  return P.each_col() % (1 / arma::sum(P, 1));

}

// [[Rcpp::export]]
arma::mat compute_probabilities_no_Gamma(const arma::mat & X, const arma::colvec & alpha, const arma::mat & Beta) {

  arma::mat P = arma::exp(compute_linear_predictor_no_Gamma(X, alpha, Beta));

  return P.each_col() % (1 / arma::sum(P, 1));

}

// [[Rcpp::export]]
arma::mat compute_conditional_probabilities(const arma::mat & Y_matrix, const arma::mat & P) {

  arma::mat C = P % Y_matrix;

  return C.each_col() % (1 / arma::sum(C, 1));

}
