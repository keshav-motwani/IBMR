#' @export
glmnet_subset = function(Y_list,
                             categories,
                             category_mappings,
                             X_list,
                             Y_list_validation = NULL,
                             category_mappings_validation = NULL,
                             X_list_validation = NULL,
                             n_lambda = 25,
                             lambda_min_ratio = 1e-4,
                             n_iter = 1e6,
                             tolerance = 1e-9) {

  names(categories) = categories
  Y_matrix_list = lapply(1:length(Y_list), function(i) create_Y_matrix(Y_list[[i]], categories, category_mappings[[i]]))

  features = colnames(X_list[[1]])

  indices_list = lapply(Y_matrix_list, function(Y) which(rowSums(Y) == 1))

  print(paste0("Keeping ", sum(sapply(indices_list, length)), " observations out of a total of ", sum(sapply(Y_list, length))))

  Y_matrix_list = mapply(Y = Y_matrix_list, indices = indices_list, FUN = function(Y, indices) Y[indices, ], SIMPLIFY = FALSE)
  Y_list = mapply(Y = Y_list, indices = indices_list, map = category_mappings, FUN = function(Y, indices, map) unlist(map[Y[indices]]), SIMPLIFY = FALSE)
  X_list = mapply(X = X_list, indices = indices_list, FUN = function(X, indices) X[indices, ], SIMPLIFY = FALSE)

  Y = unlist(Y_list)
  X = do.call(rbind, X_list)

  fit = glmnet::glmnet(x = X, y = Y, family = "multinomial", type.multinomial = "grouped", nlambda = n_lambda, lambda.min.ratio = lambda_min_ratio, trace.it = 1, thresh = tolerance, maxit = n_iter)

  model_fits = vector("list", n_lambda)

  for (i in 1:length(fit$lambda)) {

    result = list(alpha = fit$a0[categories, i],
                  Beta = do.call(cbind, lapply(fit$beta[categories], function(x) x[, i])),
                  lambda_index = i)

    names(result$alpha) = categories
    colnames(result$Beta) = categories
    rownames(result$Beta) = features

    model_fits[[i]] = result

  }

  lambda_seq = fit$lambda

  fit = list(model_fits = model_fits,
             n_lambda = n_lambda,
             categories = categories,
             lambda_sequence = lambda_seq)

  if (!is.null(Y_list_validation)) {

    validation_negative_log_likelihood = compute_tuning_performance_no_Gamma(fit, Y_list_validation, category_mappings_validation, X_list_validation)
    fit$validation_negative_log_likelihood = validation_negative_log_likelihood

    best_tuning_parameters = which_min(validation_negative_log_likelihood)[1]

    fit$best_tuning_parameters = best_tuning_parameters
    fit$best_model = fit$model_fits[[best_tuning_parameters]]

  }

  class(fit) = "glmnet_subset"

  return(fit)

}

#' @export
glmnet_relabel = function(Y_list,
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
                          n_iter = 1e6,
                          tolerance = 1e-9) {

  fit_subset = glmnet_subset(Y_list, categories, category_mappings, X_list, Y_list_validation, category_mappings_validation, X_list_validation, n_rho, rho_min_ratio, n_iter, tolerance)

  X = do.call(rbind, X_list)

  features = colnames(X_list[[1]])

  probabilities = predict_probabilities(fit_subset$best_model, X_list)
  conditional_probabilities = predict_conditional_probabilities(probabilities, Y_list, category_mappings)
  conditional_categories = predict_categories(conditional_probabilities)
  Y = unlist(conditional_categories)

  fit = glmnet::glmnet(x = X, y = Y, family = "multinomial", type.multinomial = "grouped", nlambda = n_lambda, lambda.min.ratio = lambda_min_ratio, trace.it = 1, thresh = tolerance, maxit = n_iter)

  model_fits = vector("list", n_lambda)

  for (i in 1:length(fit$lambda)) {

    result = list(alpha = fit$a0[categories, i],
                  Beta = do.call(cbind, lapply(fit$beta[categories], function(x) x[, i])),
                  lambda_index = i,
                  rho_index = fit_subset$best_tuning_parameters)

    names(result$alpha) = categories
    colnames(result$Beta) = categories
    rownames(result$Beta) = features

    model_fits[[i]] = result

  }

  fit = list(model_fits = model_fits,
             categories = categories,
             n_rho = n_rho,
             n_lambda = n_lambda)

  if (!is.null(Y_list_validation)) {

    validation_negative_log_likelihood = compute_tuning_performance_no_Gamma(fit, Y_list_validation, category_mappings_validation, X_list_validation)
    fit$validation_negative_log_likelihood = validation_negative_log_likelihood

    best_tuning_parameters = which_min(validation_negative_log_likelihood)[1]

    fit$best_tuning_parameters = best_tuning_parameters
    fit$best_model = fit$model_fits[[best_tuning_parameters]]

  }

  class(fit) = "glmnet_relabel"

  return(fit)

}

#' @export
predict_probabilities_glmnet_split = function(fit, X_list) {

  probabilities = vector("list", length(X_list))

  for (i in 1:length(X_list)) {

    X = X_list[[i]]

    final_P = matrix(0, nrow = nrow(X), ncol = length(fit$categories))

    for (k in 1:length(fit$dataset_fits)) {

      P_fine = matrix(0, nrow = nrow(X), ncol = length(fit$categories))
      colnames(P_fine) = fit$categories

      model = fit$dataset_fits[[k]]
      alpha = model$alpha
      Beta = model$Beta
      category_mapping = model$category_mapping

      P_coarse = compute_probabilities_no_Gamma(X, alpha, Beta)
      colnames(P_coarse) = names(category_mapping)

      for (coarse_category in names(category_mapping)) {

        P_fine[, category_mapping[[coarse_category]]] = P_coarse[, coarse_category] / length(category_mapping[[coarse_category]])

      }

      final_P = final_P + P_fine

    }

    final_P = final_P / length(fit$dataset_fits)

    probabilities[[i]] = final_P

  }

  return(probabilities)

}
