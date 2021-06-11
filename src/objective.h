#ifndef OBJECTIVE
#define OBJECTIVE

#include <RcppArmadillo.h>

using namespace Rcpp;

double compute_negative_log_likelihood(const List & Y_matrix_list, const List & X_list, const List & Z_list, const arma::colvec & alpha, const arma::mat & Beta, const std::vector<arma::mat> & Gamma_list, int N);
double group_lasso_penalty(const arma::mat & Beta, double lambda);
double l2_penalty(const List & Gamma_list, double rho);
double compute_objective_function(const List & Y_matrix_list, const List & X_list, const List & Z_list, const arma::colvec & alpha, const arma::mat & Beta, const std::vector<arma::mat> & Gamma_list, double lambda, double rho, int N);
double compute_negative_log_likelihood_1(const arma::mat & Y, const arma::mat & X, const arma::mat & Z, const arma::colvec & alpha, const arma::mat & Beta, const arma::mat & Gamma, int N);

#endif
