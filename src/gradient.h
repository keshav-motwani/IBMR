#ifndef GRADIENT
#define GRADIENT

#include <RcppArmadillo.h>

using namespace Rcpp;

arma::mat compute_gradient_Beta(const List & Y_matrix_list, const List & X_list, const List & Z_list, const arma::colvec & alpha, const arma::mat & Beta, const std::vector<arma::mat> & Gamma_list, int N);
arma::colvec compute_gradient_alpha(const List & Y_matrix_list, const List & X_list, const List & Z_list, const arma::colvec & alpha, const arma::mat & Beta, const std::vector<arma::mat> & Gamma_list, int N);

#endif
