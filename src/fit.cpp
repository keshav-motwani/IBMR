#include "update.h"
#include "objective.h"
#include "step_size.h"

//' @export
// [[Rcpp::export]]
List fit_Gamma_Newton(const List & Y_matrix_list, const List & X_list, const List & Z_list, double rho, int n_iter, double tolerance, arma::field<arma::mat> Gamma_list_old) {

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

  arma::field<arma::mat> Gamma_list_new = Gamma_list_old;

  for (int i = 0; i < n_iter; i++) {

    Gamma_list_new = update_Gamma_list_Newton(Y_matrix_list, X_list, Z_list, alpha, Beta, Gamma_list_old, rho, N);

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
List fit_Gamma(const List & Y_matrix_list, const List & X_list, const List & Z_list, double rho, int n_iter, double tolerance, arma::field<arma::mat> Gamma_list_old) {

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

  arma::colvec start_step_size_Gamma = compute_min_step_size_Gamma(Y_matrix_list, Z_list, rho, N) * 100000;
  arma::colvec end_step_size_Gamma(K);

  arma::colvec alpha(q);
  alpha.zeros();

  arma::mat Beta(p, q);
  Beta.zeros();

  arma::field<arma::mat> Gamma_list_new = Gamma_list_old;

  for (int i = 0; i < n_iter; i++) {

    std::tie(Gamma_list_new, end_step_size_Gamma) = update_Gamma_list(Y_matrix_list, X_list, Z_list, alpha, Beta, Gamma_list_old, rho, N, start_step_size_Gamma);

    objective(i) = compute_objective_function(Y_matrix_list, X_list, Z_list, alpha, Beta, Gamma_list_new, 0, rho, N);

    // Rcout << i << "   " << objective(i) << "\n";

    if (i > 1 && std::abs((objective(i - 1) - objective(i)) / objective(i - 1)) < tolerance) {
      break;
    }

    Gamma_list_old = Gamma_list_new;
    start_step_size_Gamma = end_step_size_Gamma * 2;

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

  double start_step_size_alpha = compute_min_step_size_alpha(Y_matrix_list) * 10000;
  double start_step_size_Beta = compute_min_step_size_Beta(Y_matrix_list, X_list, N) * 100000;
  double end_step_size_alpha;
  double end_step_size_Beta;

  arma::field<arma::mat> Gamma_list(K);
  for (R_xlen_t i = 0; i < K; i++) {
    Gamma_list(i) = arma::mat(r, q, arma::fill::zeros);
  }

  arma::colvec alpha_new = alpha_old;
  arma::mat Beta_new = Beta_old;

  for (int i = 0; i < n_iter; i++) {

    std::tie(alpha_new, end_step_size_alpha) = update_alpha(Y_matrix_list, X_list, Z_list, alpha_old, Beta_old, Gamma_list, N, start_step_size_alpha);
    if (i == 1) Rcout << "alpha " << start_step_size_alpha / end_step_size_alpha << "\n";

    // objective(2 * i) = compute_objective_function(Y_matrix_list, X_list, Z_list, alpha_new, Beta_old, Gamma_list, lambda, 0, N);

    std::tie(Beta_new, end_step_size_Beta) = update_Beta(Y_matrix_list, X_list, Z_list, alpha_new, Beta_old, Gamma_list, lambda, N, start_step_size_Beta);
    if (i == 1) Rcout << "Beta " << start_step_size_Beta / end_step_size_Beta << "\n";

    objective(2 * i + 1) = compute_objective_function(Y_matrix_list, X_list, Z_list, alpha_new, Beta_new, Gamma_list, lambda, 0, N);

    // Rcout << i << "   " << objective(2 * i) << "\n";
    // Rcout << i << "   " << objective(2 * i + 1) << "\n";

    if (i > 1 && std::abs((objective(2 * (i - 1) + 1) - objective(2 * i + 1)) / objective(2 * (i - 1) + 1)) < tolerance) {
      break;
    }

    alpha_old = alpha_new;
    Beta_old = Beta_new;
    start_step_size_alpha = end_step_size_alpha * 2;
    start_step_size_Beta = end_step_size_Beta * 2;

  }

  return Rcpp::List::create(Rcpp::Named("alpha") = alpha_new,
                            Rcpp::Named("Beta") = Beta_new,
                            Rcpp::Named("objective") = objective);

}

//' @export
// [[Rcpp::export]]
List fit_alpha(const List & Y_matrix_list, const List & X_list, const List & Z_list, int n_iter, double tolerance, arma::colvec alpha_old) {

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

  double start_step_size_alpha = compute_min_step_size_alpha(Y_matrix_list) * 10000;
  double end_step_size_alpha;

  arma::mat Beta(p, q);
  Beta.zeros();

  arma::field<arma::mat> Gamma_list(K);
  for (R_xlen_t i = 0; i < K; i++) {
    Gamma_list(i) = arma::mat(r, q, arma::fill::zeros);
  }

  arma::colvec alpha_new = alpha_old;

  for (int i = 0; i < n_iter; i++) {

    std::tie(alpha_new, end_step_size_alpha) = update_alpha(Y_matrix_list, X_list, Z_list, alpha_old, Beta, Gamma_list, N, start_step_size_alpha);
    if (i == 0) Rcout << "alpha " << start_step_size_alpha / end_step_size_alpha << "\n";

    objective(i) = compute_objective_function(Y_matrix_list, X_list, Z_list, alpha_new, Beta, Gamma_list, 0, 0, N);

    // Rcout << i << "   " << objective(i) << "\n";

    if (i > 1 && std::abs((objective(i - 1) - objective(i)) / objective(i - 1)) < tolerance) {
      break;
    }

    alpha_old = alpha_new;
    start_step_size_alpha = end_step_size_alpha * 2;

  }

  return Rcpp::List::create(Rcpp::Named("alpha") = alpha_new,
                            Rcpp::Named("objective") = objective);

}

//' @export
// [[Rcpp::export]]
List fit_alpha_Gamma_Newton(const List & Y_matrix_list, const List & X_list, const List & Z_list, double rho, int n_iter, double tolerance, arma::colvec alpha_old, arma::field<arma::mat> Gamma_list_old) {

  arma::colvec objective(2 * n_iter);
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

  double start_step_size_alpha = compute_min_step_size_alpha(Y_matrix_list) * 10000;
  double end_step_size_alpha;

  arma::mat Beta(p, q);
  Beta.zeros();

  arma::colvec alpha_new = alpha_old;
  arma::field<arma::mat> Gamma_list_new = Gamma_list_old;

  for (int i = 0; i < n_iter; i++) {

    std::tie(alpha_new, end_step_size_alpha) = update_alpha(Y_matrix_list, X_list, Z_list, alpha_old, Beta, Gamma_list_old, N, start_step_size_alpha);

    // objective(2 * i) = compute_objective_function(Y_matrix_list, X_list, Z_list, alpha_new, Beta_old, Gamma_list, lambda, 0, N);

    Gamma_list_new = update_Gamma_list_Newton(Y_matrix_list, X_list, Z_list, alpha_new, Beta, Gamma_list_old, rho, N);

    objective(2 * i + 1) = compute_objective_function(Y_matrix_list, X_list, Z_list, alpha_new, Beta, Gamma_list_new, 0, rho, N);

    // Rcout << i << "   " << objective(2 * i) << "\n";
    // Rcout << i << "   " << objective(2 * i + 1) << "\n";

    if (i > 1 && std::abs((objective(2 * (i - 1) + 1) - objective(2 * i + 1)) / objective(2 * (i - 1) + 1)) < tolerance) {
      break;
    }

    alpha_old = alpha_new;
    Gamma_list_old = Gamma_list_new;
    start_step_size_alpha = end_step_size_alpha * 2;

  }

  return Rcpp::List::create(Rcpp::Named("alpha") = alpha_new,
                            Rcpp::Named("Gamma_list") = Gamma_list_new,
                            Rcpp::Named("objective") = objective);

}

//' @export
// [[Rcpp::export]]
List fit_alpha_Gamma(const List & Y_matrix_list, const List & X_list, const List & Z_list, double rho, int n_iter, double tolerance, arma::colvec alpha_old, arma::field<arma::mat> Gamma_list_old) {

  arma::colvec objective(2 * n_iter);
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

  double start_step_size_alpha = compute_min_step_size_alpha(Y_matrix_list) * 10000;
  arma::colvec start_step_size_Gamma = compute_min_step_size_Gamma(Y_matrix_list, Z_list, rho, N) * 100000;
  double end_step_size_alpha;
  arma::colvec end_step_size_Gamma(K);

  arma::mat Beta(p, q);
  Beta.zeros();

  arma::colvec alpha_new = alpha_old;
  arma::field<arma::mat> Gamma_list_new = Gamma_list_old;

  for (int i = 0; i < n_iter; i++) {

    std::tie(alpha_new, end_step_size_alpha) = update_alpha(Y_matrix_list, X_list, Z_list, alpha_old, Beta, Gamma_list_old, N, start_step_size_alpha);
    if (i == 0) Rcout << "alpha " << start_step_size_alpha / end_step_size_alpha << "\n";

    // objective(2 * i) = compute_objective_function(Y_matrix_list, X_list, Z_list, alpha_new, Beta_old, Gamma_list, lambda, 0, N);

    std::tie(Gamma_list_new, end_step_size_Gamma) = update_Gamma_list(Y_matrix_list, X_list, Z_list, alpha_new, Beta, Gamma_list_old, rho, N, start_step_size_Gamma);
    if (i == 0) Rcout << "Gamma " << start_step_size_Gamma / end_step_size_Gamma << "\n";

    objective(2 * i + 1) = compute_objective_function(Y_matrix_list, X_list, Z_list, alpha_new, Beta, Gamma_list_new, 0, rho, N);

    // Rcout << i << "   " << objective(2 * i) << "\n";
    // Rcout << i << "   " << objective(2 * i + 1) << "\n";

    if (i > 1 && std::abs((objective(2 * (i - 1) + 1) - objective(2 * i + 1)) / objective(2 * (i - 1) + 1)) < tolerance) {
      break;
    }

    alpha_old = alpha_new;
    Gamma_list_old = Gamma_list_new;
    start_step_size_alpha = end_step_size_alpha * 2;
    start_step_size_Gamma = end_step_size_Gamma * 2;

  }

  return Rcpp::List::create(Rcpp::Named("alpha") = alpha_new,
                            Rcpp::Named("Gamma_list") = Gamma_list_new,
                            Rcpp::Named("objective") = objective);

}

//' @export
// [[Rcpp::export]]
List fit_alpha_Beta_Gamma_Newton(const List & Y_matrix_list, const List & X_list, const List & Z_list, double lambda, double rho, int n_iter, double tolerance, arma::colvec alpha_old, arma::mat Beta_old, arma::field<arma::mat> Gamma_list_old) {

  arma::colvec objective(3 * n_iter);
  objective.zeros();

  R_xlen_t K = Y_matrix_list.size();

  int N = 0;
  for (R_xlen_t i = 0; i < K; i++) {
    NumericMatrix Y = Y_matrix_list[i];
    N += Y.nrow();
  }

  double start_step_size_alpha = compute_min_step_size_alpha(Y_matrix_list) * 10000;
  double start_step_size_Beta = compute_min_step_size_Beta(Y_matrix_list, X_list, N) * 100000;
  double end_step_size_alpha;
  double end_step_size_Beta;

  arma::colvec alpha_new = alpha_old;
  arma::mat Beta_new = Beta_old;
  arma::field<arma::mat> Gamma_list_new = Gamma_list_old;

  for (int i = 0; i < n_iter; i++) {

    std::tie(alpha_new, end_step_size_alpha) = update_alpha(Y_matrix_list, X_list, Z_list, alpha_old, Beta_old, Gamma_list_old, N, start_step_size_alpha);

    // objective(3 * i) = compute_objective_function(Y_matrix_list, X_list, Z_list, alpha_new, Beta_old, Gamma_list_old, lambda, rho, N);

    std::tie(Beta_new, end_step_size_Beta) = update_Beta(Y_matrix_list, X_list, Z_list, alpha_new, Beta_old, Gamma_list_old, lambda, N, start_step_size_Beta);

    // objective(3 * i + 1) = compute_objective_function(Y_matrix_list, X_list, Z_list, alpha_new, Beta_new, Gamma_list_old, lambda, rho, N);

    Gamma_list_new = update_Gamma_list_Newton(Y_matrix_list, X_list, Z_list, alpha_new, Beta_new, Gamma_list_old, rho, N);

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
    start_step_size_alpha = end_step_size_alpha * 2;
    start_step_size_Beta = end_step_size_Beta * 2;

  }

  return Rcpp::List::create(Rcpp::Named("alpha") = alpha_new,
                            Rcpp::Named("Beta") = Beta_new,
                            Rcpp::Named("Gamma_list") = Gamma_list_new,
                            Rcpp::Named("objective") = objective);

}

//' @export
// [[Rcpp::export]]
List fit_alpha_Beta_Gamma(const List & Y_matrix_list, const List & X_list, const List & Z_list, double lambda, double rho, int n_iter, double tolerance, arma::colvec alpha_old, arma::mat Beta_old, arma::field<arma::mat> Gamma_list_old) {

  arma::colvec objective(3 * n_iter);
  objective.zeros();

  R_xlen_t K = Y_matrix_list.size();

  int N = 0;
  for (R_xlen_t i = 0; i < K; i++) {
    NumericMatrix Y = Y_matrix_list[i];
    N += Y.nrow();
  }

  double start_step_size_alpha = compute_min_step_size_alpha(Y_matrix_list) * 10000;
  double start_step_size_Beta = compute_min_step_size_Beta(Y_matrix_list, X_list, N) * 100000;
  arma::colvec start_step_size_Gamma = compute_min_step_size_Gamma(Y_matrix_list, Z_list, rho, N) * 100000;
  double end_step_size_alpha;
  double end_step_size_Beta;
  arma::colvec end_step_size_Gamma(K);

  arma::colvec alpha_new = alpha_old;
  arma::mat Beta_new = Beta_old;
  arma::field<arma::mat> Gamma_list_new = Gamma_list_old;

  for (int i = 0; i < n_iter; i++) {

    std::tie(alpha_new, end_step_size_alpha) = update_alpha(Y_matrix_list, X_list, Z_list, alpha_old, Beta_old, Gamma_list_old, N, start_step_size_alpha);
    if (i == 0) Rcout << "alpha " << start_step_size_alpha / end_step_size_alpha << "\n";

    // objective(3 * i) = compute_objective_function(Y_matrix_list, X_list, Z_list, alpha_new, Beta_old, Gamma_list_old, lambda, rho, N);

    std::tie(Beta_new, end_step_size_Beta) = update_Beta(Y_matrix_list, X_list, Z_list, alpha_new, Beta_old, Gamma_list_old, lambda, N, start_step_size_Beta);
    if (i == 0) Rcout << "Beta " << start_step_size_Beta / end_step_size_Beta << "\n";

    // objective(3 * i + 1) = compute_objective_function(Y_matrix_list, X_list, Z_list, alpha_new, Beta_new, Gamma_list_old, lambda, rho, N);

    std::tie(Gamma_list_new, end_step_size_Gamma) = update_Gamma_list(Y_matrix_list, X_list, Z_list, alpha_new, Beta_new, Gamma_list_old, rho, N, start_step_size_Gamma);
    if (i == 0) Rcout << "Gamma " << start_step_size_Gamma / end_step_size_Gamma << "\n";

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
    start_step_size_alpha = end_step_size_alpha * 2;
    start_step_size_Beta = end_step_size_Beta * 2;
    start_step_size_Gamma = end_step_size_Gamma * 2;

  }

  return Rcpp::List::create(Rcpp::Named("alpha") = alpha_new,
                            Rcpp::Named("Beta") = Beta_new,
                            Rcpp::Named("Gamma_list") = Gamma_list_new,
                            Rcpp::Named("objective") = objective);

}
