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
                             n_iter = 1e5,
                             tolerance = 1e-7) {

  names(categories) = categories
  Y_matrix_list = lapply(1:length(Y_list), function(i) create_Y_matrix(Y_list[[i]], categories, as.list(categories)))

  features = colnames(X_list[[1]])

  indices_list = lapply(Y_matrix_list, function(Y) which(rowSums(Y) == 1))
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
             n_lambda = length(fit$lambda),
             categories = categories,
             lambda_sequence = lambda_seq)

  if (!is.null(Y_list_validation)) {

    validation_negative_log_likelihood = compute_tuning_performance_no_Gamma(fit, Y_list_validation, category_mappings_validation, X_list_validation)
    fit$validation_negative_log_likelihood = validation_negative_log_likelihood

    best_tuning_parameters = which_min(validation_negative_log_likelihood)

    fit$best_tuning_parameters = best_tuning_parameters
    fit$best_model = fit$model_fits[[best_tuning_parameters]]

  }

  class(fit) = "glmnet_subset"

  return(fit)

}

#' @export
glmnet_split = function(Y_list,
                            categories,
                            category_mappings,
                            X_list,
                            Y_list_validation = NULL,
                            category_mappings_validation = NULL,
                            X_list_validation = NULL,
                            n_lambda = 25,
                            lambda_min_ratio = 1e-4,
                            n_iter = 1e5,
                            tolerance = 1e-7) {

  Y_list = c(Y_list, Y_list_validation)
  category_mappings = c(category_mappings, category_mappings_validation)
  X_list = c(X_list, X_list_validation)

  features = colnames(X_list[[1]])

  dataset_fits = vector("list", length(Y_list))

  for (k in 1:length(Y_list)) {

    fit = glmnet::cv.glmnet(x = X_list[[k]], y = Y_list[[k]], family = "multinomial", type.multinomial = "grouped", nlambda = n_lambda, lambda.min.ratio = lambda_min_ratio, trace.it = 1, thresh = tolerance, maxit = n_iter)

    lambda_index = which(fit$lambda == fit$lambda.min)

    result = list(alpha = fit$glmnet.fit$a0[names(category_mappings[[k]]), lambda_index],
                  Beta = do.call(cbind, lapply(fit$glmnet.fit$beta[names(category_mappings[[k]])], function(x) x[, lambda_index])),
                  lambda_index = lambda_index,
                  category_mapping = category_mappings[[k]])

    names(result$alpha) = names(category_mappings[[k]])
    colnames(result$Beta) = names(category_mappings[[k]])
    rownames(result$Beta) = features

    dataset_fits[[k]] = result

  }

  fit = list(dataset_fits = dataset_fits,
             categories = categories)

  class(fit) = "glmnet_split"

  return(fit)

}

#' @export
glmnet_relabel = function(Y_list,
                              categories,
                              category_mappings,
                              X_list,
                              Y_list_validation = NULL,
                              category_mappings_validation = NULL,
                              X_list_validation = NULL,
                              n_rho = 25,
                              rho_min_ratio = 1e-4,
                              n_lambda = 25,
                              lambda_min_ratio = 1e-4,
                              n_iter = 1e5,
                              tolerance = 1e-7) {

  fit_subset = glmnet_subset(Y_list, categories, category_mappings, X_list, Y_list_validation, category_mappings_validation, X_list_validation, n_rho, rho_min_ratio, n_iter, tolerance)

  model_fits = vector("list", n_rho)

  X = do.call(rbind, X_list)

  features = colnames(X_list[[1]])

  for (r in 1:fit_subset$n_lambda) {

    probabilities = predict_probabilities(fit_subset$model_fits[[r]], X_list)
    conditional_probabilities = predict_conditional_probabilities(probabilities, Y_list, category_mappings)
    conditional_categories = predict_categories(conditional_probabilities)
    Y = unlist(conditional_categories)

    fit = glmnet::glmnet(x = X, y = Y, family = "multinomial", type.multinomial = "grouped", nlambda = n_lambda, lambda.min.ratio = lambda_min_ratio, trace.it = 1, thresh = tolerance, maxit = n_iter)

    model_fits_fixed_rho = vector("list", n_lambda)

    for (i in 1:length(fit$lambda)) {

      result = list(alpha = fit$a0[categories, i],
                    Beta = do.call(cbind, lapply(fit$beta[categories], function(x) x[, i])),
                    lambda_index = i,
                    rho_index = r)

      names(result$alpha) = categories
      colnames(result$Beta) = categories
      rownames(result$Beta) = features

      model_fits_fixed_rho[[i]] = result

    }

    model_fits[[r]] = model_fits_fixed_rho

  }

  fit = list(model_fits = model_fits,
             categories = categories,
             n_rho = n_rho,
             n_lambda = n_lambda)

  if (!is.null(Y_list_validation)) {

    validation_negative_log_likelihood = compute_tuning_performance(fit, Y_list_validation, category_mappings_validation, X_list_validation)
    fit$validation_negative_log_likelihood = validation_negative_log_likelihood

    best_tuning_parameters = which_min(validation_negative_log_likelihood)

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
