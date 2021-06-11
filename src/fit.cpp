#include "update.h"
#include "objective.h"
#include "step_size.h"

//' @export
// [[Rcpp::export]]
List fit_Gamma_fast(const List & Y_matrix_list, const List & X_list, const List & Z_list, double rho, int n_iter, double tolerance, std::vector<arma::mat> Gamma_list_old) {

  arma::colvec objective(n_iter);
  objective.zeros();

  R_xlen_t K = Y_matrix_list.size();

  int N = 0;
  for (R_xlen_t i = 0; i < K; i++) {
    NumericMatrix Y = Y_matrix_list[i];
    N += Y.nrow();
  }

  NumericMatrix X0 = X_list[0];
  NumericMatrix Y0 = Y_matrix_list[0];
  int p = X0.ncol();
  int q = Y0.ncol();

  arma::colvec alpha(q);
  alpha.zeros();

  arma::mat Beta(p, q);
  Beta.zeros();

  std::vector<arma::mat> Gamma_list_new = Gamma_list_old;

  for (int i = 0; i < n_iter; i++) {

    Gamma_list_new = update_Gamma_list_fast(Y_matrix_list, X_list, Z_list, alpha, Beta, Gamma_list_old, rho, N);

    objective(i) = compute_objective_function(Y_matrix_list, X_list, Z_list, alpha, Beta, Gamma_list_new, 0, rho, N);

    // Rcout << i << "   " << objective(i) << "\n";

    if (i > 1 && std::abs((objective(i - 1) - objective(i)) / objective(i - 1)) < tolerance) {
      break;
    }

    Gamma_list_old = Gamma_list_new;

  }

  return Rcpp::List::create(Rcpp::Named("Gamma_list") = Gamma_list_new,
                            Rcpp::Named("objective") = objective);
}

//' @export
// [[Rcpp::export]]
List fit_Gamma(const List & Y_matrix_list, const List & X_list, const List & Z_list, double rho, int n_iter, double tolerance, std::vector<arma::mat> Gamma_list_old) {

  arma::colvec objective(n_iter);
  objective.zeros();

  R_xlen_t K = Y_matrix_list.size();

  int N = 0;
  for (R_xlen_t i = 0; i < K; i++) {
    NumericMatrix Y = Y_matrix_list[i];
    N += Y.nrow();
  }

  NumericMatrix X0 = X_list[0];
  NumericMatrix Y0 = Y_matrix_list[0];
  int p = X0.ncol();
  int q = Y0.ncol();

  arma::colvec min_step_size_Gamma = compute_min_step_size_Gamma(Y_matrix_list, Z_list, rho, N);

  arma::colvec alpha(q);
  alpha.zeros();

  arma::mat Beta(p, q);
  Beta.zeros();

  std::vector<arma::mat> Gamma_list_new = Gamma_list_old;

  for (int i = 0; i < n_iter; i++) {

    Gamma_list_new = update_Gamma_list(Y_matrix_list, X_list, Z_list, alpha, Beta, Gamma_list_old, rho, N, min_step_size_Gamma);

    objective(i) = compute_objective_function(Y_matrix_list, X_list, Z_list, alpha, Beta, Gamma_list_new, 0, rho, N);

    // Rcout << i << "   " << objective(i) << "\n";

    if (i > 1 && std::abs((objective(i - 1) - objective(i)) / objective(i - 1)) < tolerance) {
      break;
    }

    Gamma_list_old = Gamma_list_new;

  }

  return Rcpp::List::create(Rcpp::Named("Gamma_list") = Gamma_list_new,
                            Rcpp::Named("objective") = objective);
}

//' @export
// [[Rcpp::export]]
List fit_alpha_Beta(const List & Y_matrix_list, const List & X_list, const List & Z_list, double lambda, int n_iter, double tolerance, arma::colvec alpha_old, arma::mat Beta_old) {

  arma::colvec objective(2 * n_iter);
  objective.zeros();

  R_xlen_t K = Y_matrix_list.size();

  int N = 0;
  for (R_xlen_t i = 0; i < K; i++) {
    NumericMatrix Y = Y_matrix_list[i];
    N += Y.nrow();
  }

  NumericMatrix Z0 = Z_list[0];
  NumericMatrix Y0 = Y_matrix_list[0];
  int r = Z0.ncol();
  int q = Y0.ncol();

  double min_step_size_alpha = compute_min_step_size_alpha(Y_matrix_list);
  double min_step_size_Beta = compute_min_step_size_Beta(Y_matrix_list, X_list, N);

  std::vector<arma::mat> Gamma_list(K);
  for (R_xlen_t i = 0; i < K; i++) {
    Gamma_list[i] = arma::mat(r, q, arma::fill::zeros);
  }

  arma::colvec alpha_new = alpha_old;
  arma::mat Beta_new = Beta_old;

  for (int i = 0; i < n_iter; i++) {

    alpha_new = update_alpha(Y_matrix_list, X_list, Z_list, alpha_old, Beta_old, Gamma_list, N, min_step_size_alpha);

    // objective(2 * i) = compute_objective_function(Y_matrix_list, X_list, Z_list, alpha_new, Beta_old, Gamma_list, lambda, 0, N);

    Beta_new = update_Beta(Y_matrix_list, X_list, Z_list, alpha_new, Beta_old, Gamma_list, lambda, N, min_step_size_Beta);

    objective(2 * i + 1) = compute_objective_function(Y_matrix_list, X_list, Z_list, alpha_new, Beta_new, Gamma_list, lambda, 0, N);

    // Rcout << i << "   " << objective(2 * i) << "\n";
    // Rcout << i << "   " << objective(2 * i + 1) << "\n";

    if (i > 1 && std::abs((objective(2 * (i - 1) + 1) - objective(2 * i + 1)) / objective(2 * (i - 1) + 1)) < tolerance) {
      break;
    }

    alpha_old = alpha_new;
    Beta_old = Beta_new;

  }

  return Rcpp::List::create(Rcpp::Named("alpha") = alpha_new,
                            Rcpp::Named("Beta") = Beta_new,
                            Rcpp::Named("objective") = objective);

}

//' @export
// [[Rcpp::export]]
List fit_alpha(const List & Y_matrix_list, const List & X_list, const List & Z_list, double lambda, int n_iter, double tolerance, arma::colvec alpha_old) {

  arma::colvec objective(n_iter);
  objective.zeros();

  R_xlen_t K = Y_matrix_list.size();

  int N = 0;
  for (R_xlen_t i = 0; i < K; i++) {
    NumericMatrix Y = Y_matrix_list[i];
    N += Y.nrow();
  }

  NumericMatrix X0 = X_list[0];
  NumericMatrix Z0 = Z_list[0];
  NumericMatrix Y0 = Y_matrix_list[0];
  int p = X0.ncol();
  int r = Z0.ncol();
  int q = Y0.ncol();

  double min_step_size_alpha = compute_min_step_size_alpha(Y_matrix_list);

  arma::mat Beta(p, q);
  Beta.zeros();

  std::vector<arma::mat> Gamma_list(K);
  for (R_xlen_t i = 0; i < K; i++) {
    Gamma_list[i] = arma::mat(r, q, arma::fill::zeros);
  }

  arma::colvec alpha_new = alpha_old;

  for (int i = 0; i < n_iter; i++) {

    alpha_new = update_alpha(Y_matrix_list, X_list, Z_list, alpha_old, Beta, Gamma_list, N, min_step_size_alpha);

    objective(i) = compute_objective_function(Y_matrix_list, X_list, Z_list, alpha_new, Beta, Gamma_list, 0, 0, N);

    // Rcout << i << "   " << objective(i) << "\n";

    if (i > 1 && std::abs((objective(i - 1) - objective(i)) / objective(i - 1)) < tolerance) {
      break;
    }

    alpha_old = alpha_new;

  }

  return Rcpp::List::create(Rcpp::Named("alpha") = alpha_new,
                            Rcpp::Named("objective") = objective);

}

//' @export
// [[Rcpp::export]]
List fit_alpha_Beta_Gamma_fast(const List & Y_matrix_list, const List & X_list, const List & Z_list, double lambda, double rho, int n_iter, double tolerance, arma::colvec alpha_old, arma::mat Beta_old, std::vector<arma::mat> Gamma_list_old) {

  arma::colvec objective(3 * n_iter);
  objective.zeros();

  R_xlen_t K = Y_matrix_list.size();

  int N = 0;
  for (R_xlen_t i = 0; i < K; i++) {
    NumericMatrix Y = Y_matrix_list[i];
    N += Y.nrow();
  }

  double min_step_size_alpha = compute_min_step_size_alpha(Y_matrix_list);
  double min_step_size_Beta = compute_min_step_size_Beta(Y_matrix_list, X_list, N);

  arma::colvec alpha_new = alpha_old;
  arma::mat Beta_new = Beta_old;
  std::vector<arma::mat> Gamma_list_new = Gamma_list_old;

  for (int i = 0; i < n_iter; i++) {

    alpha_new = update_alpha(Y_matrix_list, X_list, Z_list, alpha_old, Beta_old, Gamma_list_old, N, min_step_size_alpha);

    // objective(3 * i) = compute_objective_function(Y_matrix_list, X_list, Z_list, alpha_new, Beta_old, Gamma_list_old, lambda, rho, N);

    Beta_new = update_Beta(Y_matrix_list, X_list, Z_list, alpha_new, Beta_old, Gamma_list_old, lambda, N, min_step_size_Beta);

    // objective(3 * i + 1) = compute_objective_function(Y_matrix_list, X_list, Z_list, alpha_new, Beta_new, Gamma_list_old, lambda, rho, N);

    Gamma_list_new = update_Gamma_list_fast(Y_matrix_list, X_list, Z_list, alpha_new, Beta_new, Gamma_list_old, rho, N);

    objective(3 * i + 2) = compute_objective_function(Y_matrix_list, X_list, Z_list, alpha_new, Beta_new, Gamma_list_new, lambda, rho, N);

    // Rcout << i << "   " << objective(3 * i) << "\n";
    // Rcout << i << "   " << objective(3 * i + 1) << "\n";
    // Rcout << i << "   " << objective(3 * i + 2) << "\n";

    if (i > 1 && std::abs((objective(3 * (i - 1) + 2) - objective(3 * i + 2)) / objective(3 * (i - 1) + 2)) < tolerance) {
      break;
    }

    alpha_old = alpha_new;
    Beta_old = Beta_new;
    Gamma_list_old = Gamma_list_new;

  }

  return Rcpp::List::create(Rcpp::Named("alpha") = alpha_new,
                            Rcpp::Named("Beta") = Beta_new,
                            Rcpp::Named("Gamma_list") = Gamma_list_new,
                            Rcpp::Named("objective") = objective);

}

//' @export
// [[Rcpp::export]]
List fit_alpha_Beta_Gamma(const List & Y_matrix_list, const List & X_list, const List & Z_list, double lambda, double rho, int n_iter, double tolerance, arma::colvec alpha_old, arma::mat Beta_old, std::vector<arma::mat> Gamma_list_old) {

  arma::colvec objective(3 * n_iter);
  objective.zeros();

  R_xlen_t K = Y_matrix_list.size();

  int N = 0;
  for (R_xlen_t i = 0; i < K; i++) {
    NumericMatrix Y = Y_matrix_list[i];
    N += Y.nrow();
  }

  double min_step_size_alpha = compute_min_step_size_alpha(Y_matrix_list);
  double min_step_size_Beta = compute_min_step_size_Beta(Y_matrix_list, X_list, N);
  arma::colvec min_step_size_Gamma = compute_min_step_size_Gamma(Y_matrix_list, Z_list, rho, N);

  arma::colvec alpha_new = alpha_old;
  arma::mat Beta_new = Beta_old;
  std::vector<arma::mat> Gamma_list_new = Gamma_list_old;

  for (int i = 0; i < n_iter; i++) {

    alpha_new = update_alpha(Y_matrix_list, X_list, Z_list, alpha_old, Beta_old, Gamma_list_old, N, min_step_size_alpha);

    // objective(3 * i) = compute_objective_function(Y_matrix_list, X_list, Z_list, alpha_new, Beta_old, Gamma_list_old, lambda, rho, N);

    Beta_new = update_Beta(Y_matrix_list, X_list, Z_list, alpha_new, Beta_old, Gamma_list_old, lambda, N, min_step_size_Beta);

    // objective(3 * i + 1) = compute_objective_function(Y_matrix_list, X_list, Z_list, alpha_new, Beta_new, Gamma_list_old, lambda, rho, N);

    Gamma_list_new = update_Gamma_list(Y_matrix_list, X_list, Z_list, alpha_new, Beta_new, Gamma_list_old, rho, N, min_step_size_Gamma);

    objective(3 * i + 2) = compute_objective_function(Y_matrix_list, X_list, Z_list, alpha_new, Beta_new, Gamma_list_new, lambda, rho, N);

    // Rcout << i << "   " << objective(3 * i) << "\n";
    // Rcout << i << "   " << objective(3 * i + 1) << "\n";
    // Rcout << i << "   " << objective(3 * i + 2) << "\n";

    if (i > 1 && std::abs((objective(3 * (i - 1) + 2) - objective(3 * i + 2)) / objective(3 * (i - 1) + 2)) < tolerance) {
      break;
    }

    alpha_old = alpha_new;
    Beta_old = Beta_new;
    Gamma_list_old = Gamma_list_new;

  }

  return Rcpp::List::create(Rcpp::Named("alpha") = alpha_new,
                            Rcpp::Named("Beta") = Beta_new,
                            Rcpp::Named("Gamma_list") = Gamma_list_new,
                            Rcpp::Named("objective") = objective);

}
