#' @export
fit_glmnet = function(Y_list,
                      categories,
                      X_list,
                      Y_list_validation = NULL,
                      X_list_validation = NULL,
                      n_lambda = 25,
                      lambda_min_ratio = 1e-4,
                      n_alpha = 1) {

  names(categories) = categories
  Y_matrix_list = lapply(1:length(Y_list), function(i) create_Y_matrix(Y_list[[i]], categories, as.list(categories)))

  X_list = standardize_X(X_list)
  X_mean = attr(X_list, "mean")
  X_sd = attr(X_list, "sd")

  features = colnames(X_list[[1]])

  Y = unlist(Y_list)
  X = do.call(rbind, X_list)

  model_fits = vector("list", n_alpha)

  alpha_seq = seq(1, 0, length.out = n_alpha)
  lambda_grid = matrix(nrow = n_alpha, ncol = n_lambda)

  for (a in 1:n_alpha) {

    print(a)

    fit = glmnet::glmnet(x = X, y = Y, family = "multinomial", type.multinomial = "grouped", nlambda = n_lambda, lambda.min.ratio = lambda_min_ratio, alpha = alpha_seq[a], standardize = FALSE, maxit = 1e7, thresh = 1e-8, trace.it = 1)

    model_fits_lambda_seq = vector("list", n_lambda)

    for (i in 1:length(fit$lambda)) {

      result = list(alpha = fit$a0[categories, i],
                    Beta = do.call(cbind, lapply(fit$beta[categories], function(x) x[, i])),
                    lambda_index = i,
                    alpha_index = a)

      result$KKT_check = check_KKT_IBMR_no_Gamma(Y_matrix_list, X_list, fit$lambda[i], result$alpha, result$Beta)
      print(result$KKT_check)

      model_fits_lambda_seq[[i]] = result

    }

    model_fits[[a]] = model_fits_lambda_seq

    lambda_grid[a, 1:length(fit$lambda)] = fit$lambda

  }

  fit = list(model_fits = model_fits,
             n_lambda = n_lambda,
             n_alpha = n_alpha,
             n_rho = n_alpha,
             categories = categories,
             lambda_grid = lambda_grid,
             alpha_sequence = alpha_seq)

  fit = adjust_fit(fit, categories, features, X_mean, X_sd)

  if (!is.null(Y_list_validation)) {

    validation_negative_log_likelihood = compute_tuning_performance(fit, Y_list_validation, categories, create_fine_category_mappings(categories, length(Y_list_validation)), X_list_validation)
    fit$validation_negative_log_likelihood = validation_negative_log_likelihood

    best_tuning_parameters = which_min(validation_negative_log_likelihood)[1, ]

    fit$best_tuning_parameters = best_tuning_parameters
    fit$best_model = fit$model_fits[[best_tuning_parameters[1]]][[best_tuning_parameters[2]]]

    best_tuning_parameters_group_lasso = which_min(validation_negative_log_likelihood[1, ])[1]

    fit$best_tuning_parameters_group_lasso = best_tuning_parameters_group_lasso
    fit$best_model_group_lasso = fit$model_fits[[1]][[best_tuning_parameters_group_lasso]]

  }

  return(fit)

}
