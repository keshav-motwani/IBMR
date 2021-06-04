#include "prox.h"

// [[Rcpp::export]]
arma::mat group_lasso_prox(const arma::mat & matrix, double lambda) {

  arma::colvec row_norms;
  row_norms = arma::sqrt(arma::sum(arma::square(matrix), 1));

  arma::mat thresholded(matrix.n_rows, matrix.n_cols);
  thresholded.zeros();

  arma::uvec greater = arma::find(row_norms > lambda);
  thresholded.rows(greater) = matrix.rows(greater) - (arma::mat(matrix.rows(greater)).each_col() % (lambda / row_norms.elem(greater)));

  return thresholded;

}
