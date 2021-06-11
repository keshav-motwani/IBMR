#ifndef UPDATE
#define UPDATE

#include <RcppArmadillo.h>

using namespace Rcpp;

arma::mat update_Beta(const List & Y_matrix_list, const List & X_list, const List & Z_list, const arma::colvec & alpha, const arma::mat & Beta_old, const std::vector<arma::mat> & Gamma_list, double lambda, int N, double min_step_size);
arma::vec update_alpha(const List & Y_matrix_list, const List & X_list, const List & Z_list, const arma::colvec & alpha_old, const arma::mat & Beta, const std::vector<arma::mat> & Gamma_list, int N, double min_step_size);
std::vector<arma::mat> update_Gamma_list_Newton(const List & Y_matrix_list, const List & X_list, const List & Z_list, const arma::colvec & alpha, const arma::mat & Beta, const std::vector<arma::mat> & Gamma_list_old, double rho, int N);
std::vector<arma::mat> update_Gamma_list(const List & Y_matrix_list, const List & X_list, const List & Z_list, const arma::colvec & alpha, const arma::mat & Beta, const std::vector<arma::mat> & Gamma_list_old, double rho, int N, arma::colvec min_step_size);

#endif
