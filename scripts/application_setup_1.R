fit_IBMR = function(data) {

  fit = IBMR(
    data$train$Y_list,
    data$train$category_mappings$categories,
    data$train$category_mappings$category_mappings,
    data$train$X_list,
    data$train$Z_list,
    data$validation$Y_list,
    data$validation$category_mappings$category_mappings,
    data$validation$X_list,
    n_cores = 20
  )

  return(fit)

}

fit_IBMR_int = function(data) {

  fit = IBMR(
    data$train$Y_list,
    data$train$category_mappings$categories,
    data$train$category_mappings$category_mappings,
    data$train$X_list,
    data$train$Z_list_int,
    data$validation$Y_list,
    data$validation$category_mappings$category_mappings,
    data$validation$X_list,
    n_cores = 20
  )

  return(fit)

}

fit_IBMR_common_Gamma = function(data) {

  fit = IBMR(
    data$train$Y_list,
    data$train$category_mappings$categories,
    data$train$category_mappings$category_mappings,
    data$train$X_list,
    data$train$Z_list,
    data$validation$Y_list,
    data$validation$category_mappings$category_mappings,
    data$validation$X_list,
    common_Gamma = TRUE,
    n_cores = 20
  )

  return(fit)

}

fit_IBMR_no_Gamma = function(data) {

  fit = IBMR_no_Gamma(
    data$train$Y_list,
    data$train$category_mappings$categories,
    data$train$category_mappings$category_mappings,
    data$train$X_list,
    data$validation$Y_list,
    data$validation$category_mappings$category_mappings,
    data$validation$X_list
  )

  return(fit)

}

fit_IBMR_ORC_FINE = function(data) {

  fit = IBMR(
    data$train$Y_list_fine,
    data$train$category_mappings_fine$categories,
    data$train$category_mappings_fine$category_mappings,
    data$train$X_list,
    data$train$Z_list,
    data$validation$Y_list_fine,
    data$validation$category_mappings_fine$category_mappings,
    data$validation$X_list,
    n_cores = 20
  )

  return(fit)

}

fit_IBMR_int_ORC_FINE = function(data) {

  fit = IBMR(
    data$train$Y_list_fine,
    data$train$category_mappings_fine$categories,
    data$train$category_mappings_fine$category_mappings,
    data$train$X_list,
    data$train$Z_list_int,
    data$validation$Y_list_fine,
    data$validation$category_mappings_fine$category_mappings,
    data$validation$X_list,
    n_cores = 20
  )

  return(fit)

}

fit_IBMR_common_Gamma_ORC_FINE = function(data) {

  fit = IBMR(
    data$train$Y_list_fine,
    data$train$category_mappings_fine$categories,
    data$train$category_mappings_fine$category_mappings,
    data$train$X_list,
    data$train$Z_list,
    data$validation$Y_list_fine,
    data$validation$category_mappings_fine$category_mappings,
    data$validation$X_list,
    common_Gamma = TRUE,
    n_cores = 20
  )

  return(fit)

}

fit_IBMR_no_Gamma_ORC_FINE = function(data) {

  fit = IBMR_no_Gamma(
    data$train$Y_list_fine,
    data$train$category_mappings_fine$categories,
    data$train$category_mappings_fine$category_mappings,
    data$train$X_list,
    data$validation$Y_list_fine,
    data$validation$category_mappings_fine$category_mappings,
    data$validation$X_list
  )

  return(fit)

}

compute_performance = function(fit, data) {

  Y_list_test = data$test$Y_list_fine
  category_mappings_test = data$test$category_mappings_fine$category_mappings
  X_list_test = data$test$X_list

  categories = fit$categories
  alpha_hat = fit$best_model$alpha
  Beta_hat = fit$best_model$Beta

  Y_matrix_list_test = lapply(1:length(Y_list_test), function(i) create_Y_matrix(Y_list_test[[i]], categories, category_mappings_test[[i]]))

  nll = compute_negative_log_likelihood_no_Gamma(Y_matrix_list_test, X_list_test, alpha_hat, Beta_hat, sum(sapply(Y_matrix_list_test, nrow)))

  error_per_dataset = sapply(1:length(Y_list_test), function(k) error(predict_categories(X_list_test[[k]], alpha_hat, Beta_hat, categories), Y_list_test[[k]]))
  mean_error = mean(error_per_dataset)

  return(list(nll = nll, mean_error = mean_error, error_per_dataset = error_per_dataset))

}

evaluate_parameters = function(data, method) {

  method_function = get(paste0("fit_", method))

  fit = method_function(data)

  performance = compute_performance(fit, data)

  return(list(method = method, fit = fit, performance = performance))

}
