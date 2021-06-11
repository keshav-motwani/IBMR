#ifndef STEP_SIZE
#define STEP_SIZE

#include <RcppArmadillo.h>

using namespace Rcpp;

double compute_min_step_size_Beta(const List & Y_matrix_list, const List & X_list, int N);
double compute_min_step_size_alpha(const List & Y_matrix_list);
arma::colvec compute_min_step_size_Gamma(const List & Y_matrix_list, const List & Z_list, double rho, int N);

#endif
