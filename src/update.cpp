#include "gradient.h"
#include "objective.h"
#include "prox.h"
#include "probabilities.h"
#include "update.h"

// [[Rcpp::export]]
arma::mat update_Beta(const List & Y_matrix_list, const List & X_list, const List & Z_list, const arma::colvec & alpha, const arma::mat & Beta_old, const std::vector<arma::mat> & Gamma_list, double lambda, int N) {

  bool line_search = true;
  double step_size = 0.1;
  double shrinkage = 0.5;

  arma::mat Beta_new;
  arma::mat gradient = compute_gradient_Beta(Y_matrix_list, X_list, Z_list, alpha, Beta_old, Gamma_list, N);

  double g_old = compute_negative_log_likelihood(Y_matrix_list, X_list, Z_list, alpha, Beta_old, Gamma_list, N);

  double g_new;

  while(line_search) {

    Beta_new = group_lasso_prox(Beta_old - step_size * gradient, step_size * lambda);
    g_new = compute_negative_log_likelihood(Y_matrix_list, X_list, Z_list, alpha, Beta_new, Gamma_list, N);
    arma::mat difference = (Beta_old - Beta_new) / step_size;

    if (g_new > g_old - step_size * arma::accu(gradient % difference) + 0.5 * step_size * arma::accu(arma::pow(difference, 2))) {
      step_size = shrinkage * step_size;
    } else {
      line_search = false;
    }

  }

  return Beta_new;

}

// [[Rcpp::export]]
arma::vec update_alpha(const List & Y_matrix_list, const List & X_list, const List & Z_list, const arma::colvec & alpha_old, const arma::mat & Beta, const std::vector<arma::mat> & Gamma_list, int N) {

  bool line_search = true;
  double step_size = 0.01;
  double shrinkage = 0.5;

  arma::colvec alpha_new;
  arma::colvec gradient = compute_gradient_alpha(Y_matrix_list, X_list, Z_list, alpha_old, Beta, Gamma_list, N);
  double g_old = compute_negative_log_likelihood(Y_matrix_list, X_list, Z_list, alpha_old, Beta, Gamma_list, N);

  while(line_search) {

    alpha_new = alpha_old - step_size * gradient;
    double g_new = compute_negative_log_likelihood(Y_matrix_list, X_list, Z_list, alpha_old, Beta, Gamma_list, N);

    if (g_new > g_old - 0.5 * step_size * arma::accu(arma::pow(gradient, 2))) {
      step_size = shrinkage * step_size;
    } else {
      line_search = false;
    }

  }

  return alpha_new;

}

// [[Rcpp::export]]
std::vector<arma::mat> update_Gamma_list(const List & Y_matrix_list, const List & X_list, const List & Z_list, const arma::colvec & alpha, const arma::mat & Beta, const std::vector<arma::mat> & Gamma_list_old, double rho, int N) {

  R_xlen_t K = Y_matrix_list.size();

  std::vector<arma::mat> Gamma_list_new(K);

  for (R_xlen_t i = 0; i < K; i++) {

    NumericMatrix Y = Y_matrix_list[i];
    NumericMatrix X = X_list[i];
    NumericMatrix Z = Z_list[i];
    const arma::mat & Gamma_old = Gamma_list_old[i];

    arma::mat Y_(Y.begin(), Y.nrow(), Y.ncol(), false);
    arma::mat X_(X.begin(), X.nrow(), X.ncol(), false);
    arma::mat Z_(Z.begin(), Z.nrow(), Z.ncol(), false);

    arma::mat Gamma_new(Gamma_old);

    for (int l = 0; l < Y_.n_cols; l++) {

      arma::mat P = compute_probabilities(X_, Z_, alpha, Beta, Gamma_new);
      arma::mat C = compute_conditional_probabilities(Y_, P);

      arma::colvec S = (P.col(l) % (1 - P.col(l))) + (C.col(l) % (C.col(l) - 1));
      arma::colvec v = Z_ * Gamma_old.col(l) - ((P.col(l) - C.col(l)) / S);

      Gamma_new.col(l) = arma::solve((((Z_.each_col() % S).t() * Z_) / N) + rho * arma::eye(Z_.n_cols, Z_.n_cols), Z_.t() * (S % v) / N);

    }

    Gamma_list_new[i] = Gamma_new;

  }

  return Gamma_list_new;

}
