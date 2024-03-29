# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' @export
fit_Gamma_Newton <- function(Y_matrix_list, X_list, Z_list, rho, n_iter, tolerance, Gamma_list_old) {
    .Call(`_IBMR_fit_Gamma_Newton`, Y_matrix_list, X_list, Z_list, rho, n_iter, tolerance, Gamma_list_old)
}

#' @export
fit_Gamma <- function(Y_matrix_list, X_list, Z_list, rho, n_iter, tolerance, Gamma_list_old) {
    .Call(`_IBMR_fit_Gamma`, Y_matrix_list, X_list, Z_list, rho, n_iter, tolerance, Gamma_list_old)
}

#' @export
fit_alpha_Beta <- function(Y_matrix_list, X_list, Z_list, lambda, n_iter, tolerance, alpha_old, Beta_old) {
    .Call(`_IBMR_fit_alpha_Beta`, Y_matrix_list, X_list, Z_list, lambda, n_iter, tolerance, alpha_old, Beta_old)
}

#' @export
fit_alpha <- function(Y_matrix_list, X_list, Z_list, n_iter, tolerance, alpha_old) {
    .Call(`_IBMR_fit_alpha`, Y_matrix_list, X_list, Z_list, n_iter, tolerance, alpha_old)
}

#' @export
fit_alpha_Gamma_Newton <- function(Y_matrix_list, X_list, Z_list, rho, n_iter, tolerance, alpha_old, Gamma_list_old) {
    .Call(`_IBMR_fit_alpha_Gamma_Newton`, Y_matrix_list, X_list, Z_list, rho, n_iter, tolerance, alpha_old, Gamma_list_old)
}

#' @export
fit_alpha_Gamma <- function(Y_matrix_list, X_list, Z_list, rho, n_iter, tolerance, alpha_old, Gamma_list_old) {
    .Call(`_IBMR_fit_alpha_Gamma`, Y_matrix_list, X_list, Z_list, rho, n_iter, tolerance, alpha_old, Gamma_list_old)
}

#' @export
fit_alpha_Beta_Gamma_Newton <- function(Y_matrix_list, X_list, Z_list, lambda, rho, n_iter, tolerance, alpha_old, Beta_old, Gamma_list_old) {
    .Call(`_IBMR_fit_alpha_Beta_Gamma_Newton`, Y_matrix_list, X_list, Z_list, lambda, rho, n_iter, tolerance, alpha_old, Beta_old, Gamma_list_old)
}

#' @export
fit_alpha_Beta_Gamma <- function(Y_matrix_list, X_list, Z_list, lambda, rho, n_iter, tolerance, alpha_old, Beta_old, Gamma_list_old) {
    .Call(`_IBMR_fit_alpha_Beta_Gamma`, Y_matrix_list, X_list, Z_list, lambda, rho, n_iter, tolerance, alpha_old, Beta_old, Gamma_list_old)
}

compute_gradient_Beta <- function(Y_matrix_list, X_list, Z_list, alpha, Beta, Gamma_list, N) {
    .Call(`_IBMR_compute_gradient_Beta`, Y_matrix_list, X_list, Z_list, alpha, Beta, Gamma_list, N)
}

compute_gradient_alpha <- function(Y_matrix_list, X_list, Z_list, alpha, Beta, Gamma_list, N) {
    .Call(`_IBMR_compute_gradient_alpha`, Y_matrix_list, X_list, Z_list, alpha, Beta, Gamma_list, N)
}

compute_gradient_Gamma <- function(Y, X, Z, alpha, Beta, Gamma, rho, N) {
    .Call(`_IBMR_compute_gradient_Gamma`, Y, X, Z, alpha, Beta, Gamma, rho, N)
}

compute_negative_log_likelihood <- function(Y_matrix_list, X_list, Z_list, alpha, Beta, Gamma_list, N) {
    .Call(`_IBMR_compute_negative_log_likelihood`, Y_matrix_list, X_list, Z_list, alpha, Beta, Gamma_list, N)
}

#' @export
compute_negative_log_likelihood_no_Gamma <- function(Y_matrix_list, X_list, alpha, Beta, N) {
    .Call(`_IBMR_compute_negative_log_likelihood_no_Gamma`, Y_matrix_list, X_list, alpha, Beta, N)
}

#' @export
compute_negative_log_likelihood_from_probabilities <- function(Y_matrix_list, P_list, N) {
    .Call(`_IBMR_compute_negative_log_likelihood_from_probabilities`, Y_matrix_list, P_list, N)
}

compute_negative_log_likelihood_1 <- function(Y, X, Z, alpha, Beta, Gamma, N) {
    .Call(`_IBMR_compute_negative_log_likelihood_1`, Y, X, Z, alpha, Beta, Gamma, N)
}

group_lasso_penalty <- function(Beta, lambda) {
    .Call(`_IBMR_group_lasso_penalty`, Beta, lambda)
}

l2_penalty <- function(Gamma_list, rho) {
    .Call(`_IBMR_l2_penalty`, Gamma_list, rho)
}

compute_objective_function <- function(Y_matrix_list, X_list, Z_list, alpha, Beta, Gamma_list, lambda, rho, N) {
    .Call(`_IBMR_compute_objective_function`, Y_matrix_list, X_list, Z_list, alpha, Beta, Gamma_list, lambda, rho, N)
}

compute_linear_predictor <- function(X, Z, alpha, Beta, Gamma) {
    .Call(`_IBMR_compute_linear_predictor`, X, Z, alpha, Beta, Gamma)
}

compute_linear_predictor_no_Gamma <- function(X, alpha, Beta) {
    .Call(`_IBMR_compute_linear_predictor_no_Gamma`, X, alpha, Beta)
}

compute_probabilities <- function(X, Z, alpha, Beta, Gamma) {
    .Call(`_IBMR_compute_probabilities`, X, Z, alpha, Beta, Gamma)
}

compute_probabilities_no_Gamma <- function(X, alpha, Beta) {
    .Call(`_IBMR_compute_probabilities_no_Gamma`, X, alpha, Beta)
}

compute_conditional_probabilities <- function(Y_matrix, P) {
    .Call(`_IBMR_compute_conditional_probabilities`, Y_matrix, P)
}

group_lasso_prox <- function(matrix, lambda) {
    .Call(`_IBMR_group_lasso_prox`, matrix, lambda)
}

ridge_prox <- function(matrix, rho) {
    .Call(`_IBMR_ridge_prox`, matrix, rho)
}

compute_min_step_size_Beta <- function(Y_matrix_list, X_list, N) {
    .Call(`_IBMR_compute_min_step_size_Beta`, Y_matrix_list, X_list, N)
}

compute_min_step_size_alpha <- function(Y_matrix_list) {
    .Call(`_IBMR_compute_min_step_size_alpha`, Y_matrix_list)
}

compute_min_step_size_Gamma <- function(Y_matrix_list, Z_list, rho, N) {
    .Call(`_IBMR_compute_min_step_size_Gamma`, Y_matrix_list, Z_list, rho, N)
}

