#' @export
IBMR = function(Y_list,
                categories,
                category_mappings,
                X_list,
                Z_list,
                Y_list_validation = NULL,
                category_mappings_validation = NULL,
                X_list_validation = NULL,
                n_rho = 20,
                rho_min_ratio = 1e-6,
                n_lambda = 25,
                lambda_min_ratio = 1e-4,
                phi = 1e-3,
                n_iter = 1e5,
                tolerance = 1e-8,
                Gamma_update = "gradient",
                common_Gamma = FALSE,
                n_cores = 1) {

  Y_matrix_list = lapply(1:length(Y_list), function(i) create_Y_matrix(Y_list[[i]], categories, category_mappings[[i]]))

  features = colnames(X_list[[1]])

  if (common_Gamma) {
    Y_matrix_list = list(do.call(rbind, Y_matrix_list))
    X_list = list(do.call(rbind, X_list))
    Z_list = list(do.call(rbind, Z_list))
  }

  X_list = standardize_X(X_list)
  X_mean = attr(X_list, "mean")
  X_sd = attr(X_list, "sd")

  Z_list = standardize_Z(Z_list)
  Z_mean = attr(Z_list, "mean")
  Z_sd = attr(Z_list, "sd")

  rho_sequence = compute_rho_sequence(Y_matrix_list, X_list, Z_list, n_rho, rho_min_ratio, phi, n_iter, tolerance)
  fitted_alpha_no_Beta_no_Gamma = rho_sequence$fitted_alpha
  rho_sequence = rho_sequence$sequence

  lambda_grid = compute_lambda_grid(Y_matrix_list, X_list, Z_list, rho_sequence, n_lambda, lambda_min_ratio, n_iter, tolerance, fitted_alpha_no_Beta_no_Gamma)
  fitted_alpha_no_Beta = lambda_grid$fitted_alpha
  fitted_Gamma_no_Beta = lambda_grid$fitted_Gamma_list
  lambda_grid = lambda_grid$grid

  if (Gamma_update == "Newton") {
    fit_function = fit_alpha_Beta_Gamma_Newton
  } else {
    fit_function = fit_alpha_Beta_Gamma
  }

  model_fits = parallel::mclapply(1:n_rho, function(r)
    fit_lambda_sequence_fixed_rho(
      Y_matrix_list,
      X_list,
      Z_list,
      lambda_grid[r,],
      rho_sequence[r],
      n_iter,
      tolerance,
      fitted_alpha_no_Beta[[r]],
      matrix(0, nrow = ncol(X_list[[1]]), ncol = length(categories)),
      fitted_Gamma_no_Beta[[r]],
      r,
      fit_function
    ), mc.cores = n_cores)

  fit = list(model_fits = model_fits,
             n_lambda = n_lambda,
             n_rho = n_rho,
             categories = categories,
             rho_sequence = rho_sequence,
             lambda_grid = lambda_grid,
             Z_mean = Z_mean,
             Z_sd = Z_sd,
             common_Gamma = common_Gamma,
             no_Gamma = FALSE)

  fit = adjust_fit(fit, categories, features, X_mean, X_sd)

  if (!is.null(Y_list_validation)) {

    validation_negative_log_likelihood = compute_tuning_performance(fit, Y_list_validation, category_mappings_validation, X_list_validation)
    fit$validation_negative_log_likelihood = validation_negative_log_likelihood

    best_tuning_parameters = which_min(validation_negative_log_likelihood)[1, ]

    fit$best_tuning_parameters = best_tuning_parameters
    fit$best_model = fit$model_fits[[best_tuning_parameters[1]]][[best_tuning_parameters[2]]]

  }

  if (common_Gamma) {
    class(fit) = "IBMR_common_Gamma"
  } else {
    class(fit) = "IBMR"
  }

  return(fit)

}

fit_lambda_sequence_fixed_rho = function(Y_matrix_list,
                                         X_list,
                                         Z_list,
                                         lambda_sequence,
                                         rho,
                                         n_iter,
                                         tolerance,
                                         alpha_old,
                                         Beta_old,
                                         Gamma_list_old,
                                         rho_index,
                                         fit_function) {

  n_lambda = length(lambda_sequence)

  model_fits_lambda_sequence = vector("list", n_lambda)

  for (l in 1:n_lambda) {

    print(c(rho_index, l))

    print(system.time({fit = fit_function(Y_matrix_list, X_list, Z_list, lambda_sequence[l], rho, n_iter, tolerance, alpha_old, Beta_old, Gamma_list_old)}))
    fit$lambda_index = l
    fit$rho_index = rho_index

    fit$objective = fit$objective[fit$objective != 0]
    if (length(fit$objective) == n_iter) {
      warning(paste0("Did not converge for ", rho_index, " value of rho and ", l, " value of lambda"))
      break
    }

    alpha_old = fit$alpha
    Beta_old = fit$Beta
    Gamma_list_old = fit$Gamma_list

    fit$KKT_check = check_KKT_IBMR(Y_matrix_list, X_list, Z_list, lambda_sequence[l], rho, fit$alpha, fit$Beta, fit$Gamma_list)
    print(fit$KKT_check)

    model_fits_lambda_sequence[[l]] = fit

  }

  return(model_fits_lambda_sequence)

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
                         n_iter = 1e5,
                         tolerance = 1e-8) {

  Y_matrix_list = lapply(1:length(Y_list), function(i) create_Y_matrix(Y_list[[i]], categories, category_mappings[[i]]))
  Z_list = lapply(Y_list, function(Y) matrix(1, nrow = length(Y), ncol = 1))

  X_list = standardize_X(X_list)
  X_mean = attr(X_list, "mean")
  X_sd = attr(X_list, "sd")

  features = colnames(X_list[[1]])

  lambda_sequence = compute_lambda_sequence_no_Gamma(Y_matrix_list, X_list, Z_list, n_lambda, lambda_min_ratio, n_iter, tolerance)
  fitted_alpha_no_Beta = lambda_sequence$fitted_alpha
  lambda_sequence = lambda_sequence$sequence

  model_fits_lambda_sequence = vector("list", n_lambda)

  alpha_old = fitted_alpha_no_Beta
  Beta_old = matrix(0, nrow = ncol(X_list[[1]]), ncol = length(categories))

  for (l in 1:n_lambda) {

    print(l)

    print(system.time({fit = fit_alpha_Beta(Y_matrix_list, X_list, Z_list, lambda_sequence[l], n_iter, tolerance, alpha_old, Beta_old)}))
    fit$lambda_index = l

    fit$objective = fit$objective[fit$objective != 0]
    if (length(fit$objective) == n_iter) {
      warning(paste0("Did not converge for ", l, " value of lambda"))
      break
    }

    alpha_old = fit$alpha
    Beta_old = fit$Beta

    fit$KKT_check = check_KKT_IBMR_no_Gamma(Y_matrix_list, X_list, lambda_sequence[l], fit$alpha, fit$Beta)
    print(fit$KKT_check)

    model_fits_lambda_sequence[[l]] = fit

  }

  fit = list(model_fits = model_fits_lambda_sequence,
             n_lambda = n_lambda,
             categories = categories,
             lambda_sequence = lambda_sequence,
             common_Gamma = FALSE,
             no_Gamma = TRUE)

  fit = adjust_fit_no_Gamma(fit, categories, features, X_mean, X_sd)

  if (!is.null(Y_list_validation)) {

    validation_negative_log_likelihood = compute_tuning_performance_no_Gamma(fit, Y_list_validation, category_mappings_validation, X_list_validation)
    fit$validation_negative_log_likelihood = validation_negative_log_likelihood

    best_tuning_parameters = which_min(validation_negative_log_likelihood)[1]

    fit$best_tuning_parameters = best_tuning_parameters
    fit$best_model = fit$model_fits[[best_tuning_parameters]]

  }

  return(fit)

}
