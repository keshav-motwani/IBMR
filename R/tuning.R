compute_rho_sequence = function(Y_matrix_list, X_list, Z_list, n_rho, rho_min_ratio, phi, n_iter, tolerance) {

  alpha = fit_alpha(Y_matrix_list, X_list, Z_list, n_iter, tolerance, rep(0, ncol(Y_matrix_list[[1]])))$alpha
  Beta = matrix(0, nrow = ncol(X_list[[1]]), ncol = ncol(Y_matrix_list[[1]]))
  Gamma_list = lapply(Z_list, function(x) matrix(0, nrow = ncol(x), ncol = ncol(Y_matrix_list[[1]])))

  N = sum(sapply(X_list, nrow))

  rho = numeric(length(Y_matrix_list))

  for (k in 1:length(Y_matrix_list)) {

    rho[k] = max(abs(compute_gradient_Gamma(Y_matrix_list[[k]], X_list[[k]], Z_list[[k]], alpha, Beta, Gamma_list[[k]], 0, N))) / phi

  }

  rho = max(rho)

  return(list(sequence = log_seq(rho, rho_min_ratio * rho, n_rho), fitted_alpha = alpha))

}

compute_lambda_grid = function(Y_matrix_list, X_list, Z_list, rho_sequence, n_lambda, lambda_min_ratio, n_iter, tolerance, alpha_old) {

  lambda_grid = matrix(nrow = length(rho_sequence), ncol = n_lambda)

  Beta = matrix(0, nrow = ncol(X_list[[1]]), ncol = ncol(Y_matrix_list[[1]]))

  N = sum(sapply(X_list, nrow))

  fitted_alpha = vector("list", length(rho_sequence))
  fitted_Gamma_list = vector("list", length(rho_sequence))

  Gamma_list_old = lapply(Z_list, function(x) matrix(0, nrow = ncol(x), ncol = ncol(Y_matrix_list[[1]])))

  for (r in 1:length(rho_sequence)) {

    print(r)

    print(system.time({alpha_Gamma = fit_alpha_Gamma(Y_matrix_list, X_list, Z_list, rho_sequence[r], n_iter, tolerance, alpha_old, Gamma_list_old)}))

    gradient = compute_gradient_Beta(Y_matrix_list, X_list, Z_list, alpha_Gamma$alpha, Beta, alpha_Gamma$Gamma_list, N)

    lambda_max = max(apply(gradient, 1, function(x) sqrt(sum(x ^ 2))))

    lambda_grid[r, ] = log_seq(lambda_max, lambda_min_ratio * lambda_max, n_lambda)

    fitted_alpha[[r]] = alpha_Gamma$alpha
    fitted_Gamma_list[[r]] = alpha_Gamma$Gamma_list

    alpha_old = alpha_Gamma$alpha
    Gamma_list_old = alpha_Gamma$Gamma_list

  }

  return(list(grid = lambda_grid, fitted_alpha = fitted_alpha, fitted_Gamma_list = fitted_Gamma_list))

}

compute_lambda_sequence_no_Gamma = function(Y_matrix_list, X_list, Z_list, n_lambda, lambda_min_ratio, n_iter, tolerance) {

  Beta = matrix(0, nrow = ncol(X_list[[1]]), ncol = ncol(Y_matrix_list[[1]]))
  Gamma_list = lapply(Z_list, function(x) matrix(0, nrow = ncol(x), ncol = ncol(Y_matrix_list[[1]])))

  N = sum(sapply(X_list, nrow))

  alpha = fit_alpha(Y_matrix_list, X_list, Z_list, n_iter, tolerance, rep(0, ncol(Y_matrix_list[[1]])))$alpha

  gradient = compute_gradient_Beta(Y_matrix_list, X_list, Z_list, alpha, Beta, Gamma_list, N)

  lambda_max = max(apply(gradient, 1, function(x) sqrt(sum(x ^ 2))))

  return(list(sequence = log_seq(lambda_max, lambda_min_ratio * lambda_max, n_lambda), fitted_alpha = alpha))

}

log_seq = function(from, to, length) {

  sequence = exp(seq(log(from), log(to), length.out = length))

  sequence[1] = from
  if (length > 1) sequence[length] = to

  sequence

}

#' @export
compute_tuning_performance = function(fit,
                                      Y_list_validation,
                                      categories,
                                      X_list_validation) {

  Y_matrix_list_validation = lapply(1:length(Y_list_validation), function(i) create_Y_matrix(Y_list_validation[[i]], fit$categories, category_mappings_validation[[i]]))

  N = sum(sapply(X_list_validation, nrow))

  nll = matrix(NA, nrow = fit$n_rho, ncol = fit$n_lambda)

  for (r in 1:fit$n_rho) {

    for (l in 1:fit$n_lambda) {

      print(c(r, l))

      model = fit$model_fits[[r]][[l]]

      if (!is.null(model)) nll[r, l] = compute_negative_log_likelihood_no_Gamma(Y_matrix_list_validation, X_list_validation, model$alpha, model$Beta, N)

    }

  }

  return(nll)

}

#' @export
compute_tuning_performance_no_Gamma = function(fit,
                                               Y_list_validation,
                                               category_mappings_validation,
                                               X_list_validation) {

  Y_matrix_list_validation = lapply(1:length(Y_list_validation), function(i) create_Y_matrix(Y_list_validation[[i]], fit$categories, category_mappings_validation[[i]]))

  N = sum(sapply(X_list_validation, nrow))

  nll = rep(NA, fit$n_lambda)

  for (l in 1:fit$n_lambda) {

    model = fit$model_fits[[l]]

    if (!is.null(model)) nll[l] = compute_negative_log_likelihood_no_Gamma(Y_matrix_list_validation, X_list_validation, model$alpha, model$Beta, N)

  }

  return(nll)

}

#' @export
which_min = function(mat) {

  which(mat == min(mat, na.rm = TRUE), arr.ind = TRUE)

}
