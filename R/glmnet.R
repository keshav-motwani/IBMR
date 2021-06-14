#' @export
fit_glmnet = function(Y_list,
                      categories,
                      X_list,
                      Y_list_validation,
                      X_list_validation,
                      n_lambda = 20,
                      lambda_min_ratio = 1e-4,
                      n_alpha = 25) {

  Y = unlist(Y_list)
  X = do.call(rbind, X_list)

  model_fits = vector("list", n_alpha)

  alpha_seq = seq(1, 0, length.out = n_alpha)
  lambda_grid = matrix(nrow = n_alpha, ncol = n_lambda)

  for (a in 1:n_alpha) {

    print(a)

    fit = glmnet::glmnet(x = X, y = Y, family = "multinomial", type.multinomial = "grouped", nlambda = n_lambda, lambda.min.ratio = lambda_min_ratio, alpha = alpha_seq[a])

    model_fits_lambda_seq = vector("list", n_lambda)

    for (i in 1:length(fit$lambda)) {

      result = list(alpha = fit$a0[categories, i],
                    Beta = do.call(cbind, lapply(fit$beta[categories], function(x) x[, i])),
                    lambda_index = i,
                    alpha_index = a)

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
