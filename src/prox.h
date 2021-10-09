#ifndef PROX
#define PROX

#include <RcppArmadillo.h>

using namespace Rcpp;

arma::mat group_lasso_prox(const arma::mat & matrix, double lambda);
arma::mat ridge_prox(const arma::mat & matrix, double rho);

#endif
