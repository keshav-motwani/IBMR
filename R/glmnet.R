#' @export
fit_glmnet = function(Y_list,
                      categories,
                      X_list,
                      Y_list_validation,
                      X_list_validation,
                      n_lambda = 25,
                      lambda_min_ratio = 1e-3,
                      n_iter = 100000,
                      tolerance = 1e-6) {

  Y = unlist(Y_list)
  X = do.call(rbind, X_list)

  fit = glmnet::glmnet(x = X, y = Y, family = "multinomial", type.multinomial = "grouped", nlambda = n_lambda, lambda.min.ratio = lambda_min_ratio)

  model_fits = list()

  for (i in 1:n_lambda) {

    result = list(alpha = fit$a0[categories, i],
                  Beta = do.call(cbind, lapply(fit$beta[categories], function(x) x[, i])),
                  lambda_index = i)

    model_fits = c(model_fits, list(result))

  }

  fit = list(model_fits = model_fits,
             n_lambda = n_lambda,
             categories = categories,
             lambda_sequence = fit$lambda)

  if (!is.null(Y_list_validation)) {

    validation_negative_log_likelihood = compute_tuning_performance_L0(fit, Y_list_validation, categories, X_list_validation)

    best_tuning_parameters = which_min(validation_negative_log_likelihood)[1]

    fit$best_tuning_parameters = best_tuning_parameters
    fit$best_model = fit$model_fits[[best_tuning_parameters]]
    fit$validation_negative_log_likelihood = validation_negative_log_likelihood

  }

  return(fit)

}

#' @export
fit_elastic_net = function(Y_list,
                           categories,
                           X_list,
                           Y_list_validation,
                           X_list_validation,
                           n_alpha = 2,
                           n_lambda = 3,
                           lambda_min_ratio = 1e-3,
                           n_iter = 100000,
                           tolerance = 1e-6) {

  Y = unlist(Y_list)
  X = do.call(rbind, X_list)

  model_fits = list()

  alpha_seq = seq(0, 1, length.out = n_alpha)

  for (a in 1:n_alpha) {

    print(a)

    fit = glmnet::glmnet(x = X, y = Y, family = "multinomial", type.multinomial = "grouped", nlambda = n_lambda, lambda.min.ratio = lambda_min_ratio, alpha = alpha_seq[a])

    model_fits_lambda_seq = list()

    for (i in 1:n_lambda) {

      result = list(alpha = fit$a0[categories, i],
                    Beta = do.call(cbind, lapply(fit$beta[categories], function(x) x[, i])),
                    lambda_index = i)

      model_fits_lambda_seq = c(model_fits_lambda_seq, list(result))

    }

    model_fits = c(model_fits, list(model_fits_lambda_seq))

  }

  fit = list(model_fits = model_fits,
             n_lambda = n_lambda,
             n_alpha = n_alpha,
             categories = categories,
             lambda_sequence = fit$lambda,
             alpha_sequence = alpha_seq)

  if (!is.null(Y_list_validation)) {

    fit$n_gamma = fit$n_alpha
    validation_negative_log_likelihood = compute_tuning_performance(fit, Y_list_validation, categories, X_list_validation)

    best_tuning_parameters = which_min(validation_negative_log_likelihood)[1, ]

    fit$best_tuning_parameters = best_tuning_parameters
    fit$best_model = fit$model_fits[[best_tuning_parameters[1]]][[best_tuning_parameters[2]]]
    fit$validation_negative_log_likelihood = validation_negative_log_likelihood

  }

  return(fit)

}
