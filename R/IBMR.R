#' @export
IBMR = function(Y_list,
                categories,
                category_mappings,
                X_list,
                Z_list,
                Y_list_validation = NULL,
                category_mappings_validation = NULL,
                X_list_validation = NULL,
                n_rho = 5,
                rho_min_ratio = 1e-4,
                n_lambda = 25,
                lambda_min_ratio = 1e-4,
                phi = 1e-3,
                n_iter = 1e4,
                tolerance = 1e-6,
                stop_solution_path = 1.01,
                Gamma_update = "gradient",
                common_Gamma = FALSE,
                n_cores = 1,
                verbose = TRUE) {

  Y_matrix_list = lapply(1:length(Y_list), function(i) create_Y_matrix(Y_list[[i]], categories, category_mappings[[i]]))
  if (!is.null(Y_list_validation)) {
    Y_matrix_list_validation = lapply(1:length(Y_list_validation), function(i) create_Y_matrix(Y_list_validation[[i]], categories, category_mappings_validation[[i]]))
    N_val = sum(sapply(X_list_validation, nrow))
  } else {
    Y_matrix_list_validation = NULL
    N_val = NULL
  }

  count = numeric(length(categories))
  names(count) = categories
  for (k in 1:length(Y_matrix_list)) {
    count = count + colSums(Y_matrix_list[[k]][rowSums(Y_matrix_list[[k]]) == 1, ])
  }
  if (verbose) print(count)
  stopifnot(all(count >= 1))

  features = colnames(X_list[[1]])

  if (common_Gamma) {
    Y_matrix_list = list(do.call(rbind, Y_matrix_list))
    X_list = list(do.call(rbind, X_list))
    Z_list = list(do.call(rbind, Z_list))
  }

  if (verbose) print("Standardizing predictors")

  X_list = standardize_X(X_list)
  X_mean = attr(X_list, "mean")
  X_sd = attr(X_list, "sd")

  Z_list = standardize_Z(Z_list)
  Z_mean = attr(Z_list, "mean")
  Z_sd = attr(Z_list, "sd")

  if (verbose) print("Computing tuning parameter sequences")

  rho_sequence = compute_rho_sequence(Y_matrix_list, X_list, Z_list, n_rho, rho_min_ratio, phi, n_iter, tolerance)
  fitted_alpha_no_Beta_no_Gamma = rho_sequence$fitted_alpha
  rho_sequence = rho_sequence$sequence

  lambda_grid = compute_lambda_grid(Y_matrix_list, X_list, Z_list, rho_sequence, n_lambda, lambda_min_ratio, n_iter, tolerance, fitted_alpha_no_Beta_no_Gamma)
  fitted_alpha_no_Beta = lambda_grid$fitted_alpha
  fitted_Gamma_no_Beta = lambda_grid$fitted_Gamma_list
  lambda_grid = lambda_grid$grid

  if (verbose) print("Fitting models")

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
      categories,
      features,
      Y_matrix_list_validation,
      X_list_validation,
      N_val,
      lambda_grid[r,],
      rho_sequence[r],
      n_iter,
      tolerance,
      stop_solution_path,
      verbose,
      fitted_alpha_no_Beta[[r]],
      matrix(0, nrow = ncol(X_list[[1]]), ncol = length(categories)),
      fitted_Gamma_no_Beta[[r]],
      r,
      X_mean,
      X_sd,
      fit_function
    ), mc.cores = n_cores)

  validation_negative_log_likelihood = t(sapply(model_fits, `[[`, 2))
  model_fits = lapply(model_fits, `[[`, 1)

  fit = list(model_fits = model_fits,
             n_lambda = n_lambda,
             n_rho = n_rho,
             categories = categories,
             rho_sequence = rho_sequence,
             lambda_grid = lambda_grid,
             X_mean = X_mean,
             X_sd = X_sd,
             Z_mean = Z_mean,
             Z_sd = Z_sd,
             common_Gamma = common_Gamma,
             no_Gamma = FALSE)

  if (!is.null(Y_list_validation)) {

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
                                         categories,
                                         features,
                                         Y_matrix_list_validation,
                                         X_list_validation,
                                         N_val,
                                         lambda_sequence,
                                         rho,
                                         n_iter,
                                         tolerance,
                                         stop_solution_path,
                                         verbose,
                                         alpha_old,
                                         Beta_old,
                                         Gamma_list_old,
                                         rho_index,
                                         X_mean,
                                         X_sd,
                                         fit_function) {

  n_lambda = length(lambda_sequence)

  validation_negative_log_likelihood = rep(NA, n_lambda)

  model_fits_lambda_sequence = vector("list", n_lambda)

  for (l in 1:n_lambda) {

    if (verbose) print(paste0("Fitting tuning parameters: ", rho_index, ", ", l))

    fit = fit_function(Y_matrix_list, X_list, Z_list, lambda_sequence[l], rho, n_iter, tolerance, alpha_old, Beta_old, Gamma_list_old)
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
    if (verbose) print(fit$KKT_check)
    fit$alpha = adjust_alpha(fit$alpha, fit$Beta, X_mean, X_sd)
    fit$Beta = adjust_Beta(fit$Beta, X_sd)
    names(fit$alpha) = categories
    colnames(fit$Beta) = categories
    rownames(fit$Beta) = features

    model_fits_lambda_sequence[[l]] = fit

    if (!is.null(Y_matrix_list_validation)) {

      validation_negative_log_likelihood[l] = compute_negative_log_likelihood_no_Gamma(Y_matrix_list_validation, X_list_validation, fit$alpha, fit$Beta, N_val)

      if ((l > n_lambda / 4) && !is.na(stop_solution_path) && validation_negative_log_likelihood[l] > stop_solution_path * min(validation_negative_log_likelihood, na.rm = TRUE)) {
        break
      }

    }

  }

  return(list(model_fits = model_fits_lambda_sequence, validation_negative_log_likelihood = validation_negative_log_likelihood))

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
                         n_iter = 1e4,
                         tolerance = 1e-6,
                         stop_solution_path = 1.01,
                         verbose = TRUE) {

  Y_matrix_list = lapply(1:length(Y_list), function(i) create_Y_matrix(Y_list[[i]], categories, category_mappings[[i]]))
  Y_matrix_list = list(do.call(rbind, Y_matrix_list))

  count = numeric(length(categories))
  names(count) = categories
  for (k in 1:length(Y_matrix_list)) {
    count = count + colSums(Y_matrix_list[[k]][rowSums(Y_matrix_list[[k]]) == 1, ])
  }
  if (verbose) print(count)
  stopifnot(all(count >= 1))

  X_list = list(do.call(rbind, X_list))
  Z_list = list(matrix(1, nrow = nrow(X_list[[1]]), ncol = 1))

  if (!is.null(Y_list_validation)) {
    Y_matrix_list_validation = lapply(1:length(Y_list_validation), function(i) create_Y_matrix(Y_list_validation[[i]], categories, category_mappings_validation[[i]]))
    N_val = sum(sapply(X_list_validation, nrow))
    validation_negative_log_likelihood = rep(NA, n_lambda)
  }

  if (verbose) print("Standardizing predictors")

  X_list = standardize_X(X_list)
  X_mean = attr(X_list, "mean")
  X_sd = attr(X_list, "sd")

  features = colnames(X_list[[1]])

  if (verbose) print("Computing tuning parameter sequences")

  lambda_sequence = compute_lambda_sequence_no_Gamma(Y_matrix_list, X_list, Z_list, n_lambda, lambda_min_ratio, n_iter, tolerance)
  fitted_alpha_no_Beta = lambda_sequence$fitted_alpha
  lambda_sequence = lambda_sequence$sequence

  if (verbose) print("Fitting models")

  model_fits_lambda_sequence = vector("list", n_lambda)

  alpha_old = fitted_alpha_no_Beta
  Beta_old = matrix(0, nrow = ncol(X_list[[1]]), ncol = length(categories))

  for (l in 1:n_lambda) {

    if (verbose) print(paste0("Fitting tuning parameter: ", l))

    fit = fit_alpha_Beta(Y_matrix_list, X_list, Z_list, lambda_sequence[l], n_iter, tolerance, alpha_old, Beta_old)
    fit$lambda_index = l

    fit$objective = fit$objective[fit$objective != 0]
    if (length(fit$objective) == n_iter) {
      warning(paste0("Did not converge for ", l, " value of lambda"))
      break
    }

    alpha_old = fit$alpha
    Beta_old = fit$Beta

    fit$KKT_check = check_KKT_IBMR_no_Gamma(Y_matrix_list, X_list, lambda_sequence[l], fit$alpha, fit$Beta)
    if (verbose) print(fit$KKT_check)

    fit$alpha = adjust_alpha(fit$alpha, fit$Beta, X_mean, X_sd)
    fit$Beta = adjust_Beta(fit$Beta, X_sd)
    names(fit$alpha) = categories
    colnames(fit$Beta) = categories
    rownames(fit$Beta) = features

    model_fits_lambda_sequence[[l]] = fit

    if (!is.null(Y_list_validation)) {

      validation_negative_log_likelihood[l] = compute_negative_log_likelihood_no_Gamma(Y_matrix_list_validation, X_list_validation, fit$alpha, fit$Beta, N_val)

      if (l > n_lambda / 4 && !is.na(stop_solution_path) && validation_negative_log_likelihood[l] > stop_solution_path * min(validation_negative_log_likelihood, na.rm = TRUE)) {
        break
      }

    }

  }

  fit = list(model_fits = model_fits_lambda_sequence,
             n_lambda = n_lambda,
             categories = categories,
             X_mean = X_mean, X_sd = X_sd,
             lambda_sequence = lambda_sequence,
             common_Gamma = FALSE,
             no_Gamma = TRUE)

  if (!is.null(Y_list_validation)) {

    fit$validation_negative_log_likelihood = validation_negative_log_likelihood

    best_tuning_parameters = which_min(validation_negative_log_likelihood)[1]

    fit$best_tuning_parameters = best_tuning_parameters
    fit$best_model = fit$model_fits[[best_tuning_parameters]]

  }

  return(fit)

}

#' @export
subset = function(Y_list,
                  categories,
                  category_mappings,
                  X_list,
                  Y_list_validation = NULL,
                  category_mappings_validation = NULL,
                  X_list_validation = NULL,
                  n_lambda = 25,
                  lambda_min_ratio = 1e-4,
                  n_iter = 1e4,
                  tolerance = 1e-6,
                  stop_solution_path = 1.01,
                  verbose = TRUE) {

  Y_matrix_list = lapply(1:length(Y_list), function(i) create_Y_matrix(Y_list[[i]], categories, category_mappings[[i]]))

  indices_list = lapply(Y_matrix_list, function(Y) which(rowSums(Y) == 1))

  print(paste0("Keeping ", sum(sapply(indices_list, length)), " observations out of a total of ", sum(sapply(Y_list, length))))

  Y_list = list(unlist(mapply(Y = Y_list, indices = indices_list, map = category_mappings, FUN = function(Y, indices, map) unlist(map[Y[indices]]), SIMPLIFY = FALSE)))
  X_list = list(do.call(rbind, mapply(X = X_list, indices = indices_list, FUN = function(X, indices) X[indices, ], SIMPLIFY = FALSE)))

  IBMR_no_Gamma(Y_list,
                categories,
                list(as.list(setNames(nm = categories))),
                X_list,
                Y_list_validation,
                category_mappings_validation,
                X_list_validation,
                n_lambda,
                lambda_min_ratio,
                n_iter,
                tolerance,
                stop_solution_path,
                verbose)

}

#' @export
IBMR_no_Gamma_subset = subset

#' @export
relabel = function(Y_list,
                   categories,
                   category_mappings,
                   X_list,
                   Y_list_validation,
                   category_mappings_validation,
                   X_list_validation,
                   n_rho = 25,
                   rho_min_ratio = 1e-4,
                   n_lambda = 25,
                   lambda_min_ratio = 1e-4,
                   n_iter = 1e4,
                   tolerance = 1e-6,
                   stop_solution_path = 1.01,
                   verbose = TRUE) {

  fit_subset = IBMR_no_Gamma_subset(Y_list, categories, category_mappings, X_list, Y_list_validation, category_mappings_validation, X_list_validation, n_rho, rho_min_ratio, n_iter, tolerance, stop_solution_path, verbose)

  probabilities = predict_probabilities(fit_subset$best_model, X_list)
  conditional_probabilities = predict_conditional_probabilities(probabilities, Y_list, category_mappings)
  Y_list = predict_categories(conditional_probabilities)

  fit = IBMR_no_Gamma(Y_list,
                categories,
                replicate(length(Y_list), as.list(setNames(nm = categories)), simplify = FALSE),
                X_list,
                Y_list_validation,
                category_mappings_validation,
                X_list_validation,
                n_lambda,
                lambda_min_ratio,
                n_iter,
                tolerance,
                stop_solution_path,
                verbose)

  attr(fit, "subset") = fit_subset

  return(fit)

}

#' @export
IBMR_no_Gamma_relabel = relabel
