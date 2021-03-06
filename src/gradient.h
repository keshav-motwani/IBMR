#ifndef GRADIENT
#define GRADIENT

#include <RcppArmadillo.h>

using namespace Rcpp;

arma::mat compute_gradient_Beta(const List & Y_matrix_list, const List & X_list, const List & Z_list, const arma::colvec & alpha, const arma::mat & Beta, const arma::field<arma::mat> & Gamma_list, int N);
arma::colvec compute_gradient_alpha(const List & Y_matrix_list, const List & X_list, const List & Z_list, const arma::colvec & alpha, const arma::mat & Beta, const arma::field<arma::mat> & Gamma_list, int N);
arma::mat compute_gradient_Gamma(const arma::mat & Y, const arma::mat & X, const arma::mat & Z, const arma::colvec & alpha, const arma::mat & Beta, const arma::mat & Gamma, double rho, int N);

#endif
