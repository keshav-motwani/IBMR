#' @export
IBMR = function(Y_list,
                categories,
                category_mappings,
                X_list,
                Z_list,
                Y_list_validation = NULL,
                category_mappings_validation = NULL,
                X_list_validation = NULL,
                n_lambda = 25,
                lambda_min_ratio = 1e-4,
                n_rho = 20,
                rho_min_ratio = 1e-8,
                phi = 1e-3,
                n_iter = 10000,
                tolerance = 1e-8,
                Gamma_update = "gradient") {

  Y_matrix_list = lapply(1:length(Y_list), function(i) create_Y_matrix(Y_list[[i]], categories, category_mappings[[i]]))

  X_list = standardize_X(X_list)
  X_mean = attr(X_list, "mean")
  X_sd = attr(X_list, "sd")

  Z_list = standardize_Z(Z_list)
  Z_mean = attr(Z_list, "mean")
  Z_sd = attr(Z_list, "sd")

  features = colnames(X_list[[1]])

  rho_sequence = compute_rho_sequence(Y_matrix_list, X_list, Z_list, n_rho, rho_min_ratio, phi, n_iter, tolerance)
  lambda_grid = compute_lambda_grid(Y_matrix_list, X_list, Z_list, rho_sequence, n_lambda, lambda_min_ratio, n_iter, tolerance)

  if (Gamma_update == "Newton") {
    fit_function = fit_alpha_Beta_Gamma_Newton
  } else {
    fit_function = fit_alpha_Beta_Gamma
  }

  model_fits = vector("list", n_rho)

  for (r in 1:n_rho) {

    model_fits_lambda_sequence = vector("list", n_lambda)

    if (r == 1) {
      alpha_old = rep(0, length(categories))
      Beta_old = matrix(0, nrow = ncol(X_list[[1]]), ncol = length(categories))
      Gamma_list_old = lapply(Z_list, function(x) matrix(0, nrow = ncol(x), ncol = length(categories)))
    } else {
      model_old = model_fits[[r - 1]][[1]]
      alpha_old = model_old$alpha
      Beta_old = model_old$Beta
      Gamma_list_old = model_old$Gamma_list
    }

    for (l in 1:n_lambda) {

      print(c(r, l))

      fit = fit_function(Y_matrix_list, X_list, Z_list, lambda_grid[r, l], rho_sequence[r], n_iter, tolerance, alpha_old, Beta_old, Gamma_list_old)
      fit$lambda_index = l
      fit$rho_index = r

      alpha_old = fit$alpha
      Beta_old = fit$Beta
      Gamma_list_old = fit$Gamma_list

      fit$KKT_check = check_KKT_IBMR(Y_matrix_list, X_list, Z_list, lambda_grid[r, l], rho_sequence[r], fit$alpha, fit$Beta, fit$Gamma_list)
      print(fit$KKT_check)

      model_fits_lambda_sequence[[l]] = fit

    }

    model_fits[[r]] = model_fits_lambda_sequence

  }

  fit = list(model_fits = model_fits,
             n_lambda = n_lambda,
             n_rho = n_rho,
             categories = categories,
             lambda_grid = lambda_grid,
             rho_sequence = rho_sequence,
             Z_mean = Z_mean,
             Z_sd = Z_sd)

  fit = adjust_fit(fit, categories, features, X_mean, X_sd)

  if (!is.null(Y_list_validation)) {

    validation_negative_log_likelihood = compute_tuning_performance(fit, Y_list_validation, categories, category_mappings_validation, X_list_validation)
    fit$validation_negative_log_likelihood = validation_negative_log_likelihood

    best_tuning_parameters = which_min(validation_negative_log_likelihood)[1, ]

    fit$best_tuning_parameters = best_tuning_parameters
    fit$best_model = fit$model_fits[[best_tuning_parameters[1]]][[best_tuning_parameters[2]]]

  }

  return(fit)

}

#' @export
IBMR_no_Gamma = function(Y_list,
                         categories,
                         category_mappings,
                         X_list,
                         Y_list_validation = NULL,
                         category_mappings_validation = NULL,
                         X_list_validation = NULL,
                         n_lambda = 25,
                         lambda_min_ratio = 1e-4,
                         n_iter = 10000,
                         tolerance = 1e-8) {

  Y_matrix_list = lapply(1:length(Y_list), function(i) create_Y_matrix(Y_list[[i]], categories, category_mappings[[i]]))
  Z_list = lapply(Y_list, function(Y) matrix(1, nrow = length(Y), ncol = 1))

  X_list = standardize_X(X_list)
  X_mean = attr(X_list, "mean")
  X_sd = attr(X_list, "sd")

  features = colnames(X_list[[1]])

  lambda_sequence = compute_lambda_sequence_no_Gamma(Y_matrix_list, X_list, Z_list, n_lambda, lambda_min_ratio, n_iter, tolerance)

  model_fits_lambda_sequence = vector("list", n_lambda)

  alpha_old = rep(0, length(categories))
  Beta_old = matrix(0, nrow = ncol(X_list[[1]]), ncol = length(categories))
  Gamma_list_old = lapply(Z_list, function(x) matrix(0, nrow = ncol(x), ncol = length(categories)))

  for (l in 1:n_lambda) {

    print(l)

    fit = fit_alpha_Beta(Y_matrix_list, X_list, Z_list, lambda_sequence[l], n_iter, tolerance, alpha_old, Beta_old)
    fit$lambda_index = l

    alpha_old = fit$alpha
    Beta_old = fit$Beta

    fit$KKT_check = check_KKT_IBMR_no_Gamma(Y_matrix_list, X_list, lambda_sequence[l], fit$alpha, fit$Beta)
    print(fit$KKT_check)

    model_fits_lambda_sequence[[l]] = fit

  }

  fit = list(model_fits = model_fits_lambda_sequence,
             n_lambda = n_lambda,
             categories = categories,
             lambda_sequence = lambda_sequence)

  fit = adjust_fit_no_Gamma(fit, categories, features, X_mean, X_sd)

  if (!is.null(Y_list_validation)) {

    validation_negative_log_likelihood = compute_tuning_performance_no_Gamma(fit, Y_list_validation, categories, category_mappings_validation, X_list_validation)
    fit$validation_negative_log_likelihood = validation_negative_log_likelihood

    best_tuning_parameters = which_min(validation_negative_log_likelihood)[1, ]

    fit$best_tuning_parameters = best_tuning_parameters
    fit$best_model = fit$model_fits[[best_tuning_parameters[2]]]

  }

  return(fit)

}
