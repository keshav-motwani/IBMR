#ifndef UPDATE
#define UPDATE

#include <RcppArmadillo.h>

using namespace Rcpp;

std::tuple<arma::mat, double> update_Beta(const List & Y_matrix_list, const List & X_list, const List & Z_list, const arma::colvec & alpha, const arma::mat & Beta_old, const arma::field<arma::mat> & Gamma_list, double lambda, int N, double start_step_size);
std::tuple<arma::vec, double> update_alpha(const List & Y_matrix_list, const List & X_list, const List & Z_list, const arma::colvec & alpha_old, const arma::mat & Beta, const arma::field<arma::mat> & Gamma_list, int N, double start_step_size);
arma::field<arma::mat> update_Gamma_list_Newton(const List & Y_matrix_list, const List & X_list, const List & Z_list, const arma::colvec & alpha, const arma::mat & Beta, const arma::field<arma::mat> & Gamma_list_old, double rho, int N);
std::tuple<arma::field<arma::mat>, arma::colvec> update_Gamma_list(const List & Y_matrix_list, const List & X_list, const List & Z_list, const arma::colvec & alpha, const arma::mat & Beta, const arma::field<arma::mat> & Gamma_list_old, double rho, int N, arma::colvec start_step_size);

#endif
