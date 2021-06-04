#ifndef PROBABILITIES
#define PROBABILITIES

#include <RcppArmadillo.h>

using namespace Rcpp;

arma::mat compute_probabilities(const arma::mat & X, const arma::mat & Z, const arma::colvec & alpha, const arma::mat & Beta, const arma::mat & Gamma);
arma::mat compute_probabilities_Gamma0(const arma::mat & X, const arma::colvec & alpha, const arma::mat & Beta);
arma::mat compute_conditional_probabilities(const arma::mat & Y_matrix, const arma::mat & P);

#endif
