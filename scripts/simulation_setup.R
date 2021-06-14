generate_simulation_data_fine_clean = function(q,
                                               N,
                                               K,
                                               p,
                                               nonzero,
                                               b,
                                               replicate) {

  set.seed(replicate, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")

  categories = paste0("C", 1:q)
  names(categories) = categories

  category_mappings = list(
    categories = categories,
    category_mappings = replicate(K, as.list(categories), simplify = FALSE),
    inverse_category_mappings = replicate(K, categories, simplify = FALSE)
  )
  category_mappings_fine = category_mappings

  alpha = simulate_alpha(category_mappings$categories)
  Beta = simulate_Beta(category_mappings$categories, p, nonzero, -b, b)

  X_star_list = simulate_X_star_list(rep(N / K, K), p)
  Y_list = simulate_Y_list(category_mappings$categories, category_mappings$inverse_category_mappings, X_star_list, alpha, Beta)
  X_list = X_star_list
  Z_list = compute_pca_for_Z_list(X_list, 50, TRUE)

  X_star_list_val = simulate_X_star_list(rep(N / K, K), p)
  Y_list_val = simulate_Y_list(category_mappings$categories, category_mappings$inverse_category_mappings, X_star_list_val, alpha, Beta)
  X_list_val = X_star_list_val

  X_star_list_test = simulate_X_star_list(N, p)
  Y_list_test = simulate_Y_list(category_mappings_fine$categories, category_mappings_fine$inverse_category_mappings[1], X_star_list_test, alpha, Beta)

  output = list(train = list(X_list = X_list,
                             X_star_list = X_star_list,
                             Z_list = Z_list,
                             Y_list = Y_list,
                             Y_list_fine = get_fine_categories(Y_list),
                             category_mappings = category_mappings,
                             category_mappings_fine = category_mappings_fine),
                validation = list(X_list = X_list_val,
                                  X_star_list = X_star_list_val,
                                  Y_list = Y_list_val,
                                  Y_list_fine = get_fine_categories(Y_list_val),
                                  category_mappings = category_mappings,
                                  category_mappings_fine = category_mappings_fine),
                test = list(X_star_list = X_star_list_test,
                            Y_list_fine = Y_list_test,
                            category_mappings_fine = category_mappings_fine),
                alpha = alpha,
                Beta = Beta)

  return(output)

}

get_fine_categories = function(Y_list) {

  lapply(Y_list, function(Y) names(Y))

}

expand_parameters = function(run_name,
                             considered_values,
                             defaults,
                             n_replicates,
                             methods) {

  parameter_list = list()

  for (name in names(considered_values)) {
    values = considered_values[[name]]
    for (replicate in 1:n_replicates) {
      for (value in values) {
        for (method in methods) {
          params = c(
            defaults,
            experiment = name,
            replicate = replicate,
            method = method,
            run = run_name
          )
          params[name] = value
          parameter_list = c(parameter_list, list(params))
        }
      }
    }
  }

  return(parameter_list)
}

evaluate_parameters = function(parameters, simulation_function) {

  method = parameters$method
  method_function = get(paste0("fit_", method))

  data = do.call(simulation_function, parameters[formalArgs(simulation_function)])

  fits = method_function(data)

  results = vector("list", length(fits))

  for (i in 1:length(fits)) {

    fit = fits[[i]]

    performance = compute_performance(data$alpha,
                                      data$Beta,
                                      fit$alpha,
                                      fit$Beta,
                                      data$test$Y_list_fine,
                                      data$test$X_star_list,
                                      data$test$category_mappings_fine$categories)
    performance$tuning_parameters = fit$tuning_parameters
    performance$validation_negative_log_likelihood = fit$validation_negative_log_likelihood

    best_case_performance = compute_best_case_performance(data$alpha,
                                                          data$Beta,
                                                          fit$all_alphas,
                                                          fit$all_Betas,
                                                          data$test$Y_list_fine,
                                                          data$test$X_star_list,
                                                          data$test$category_mappings_fine$categories)

    if (length(fits) > 1) parameters$method = names(fits)[i]

    results[[i]] = list(parameters = parameters, performance = performance, best_case_performance = best_case_performance)

  }

  return(results)

}

compute_performance = function(alpha, Beta, alpha_hat, Beta_hat, Y_list_test, X_list_test, categories) {

  Beta_SSE = Beta_SSE(Beta_hat, Beta)
  Beta_FPR = Beta_FPR(Beta_hat, Beta)
  Beta_TPR = Beta_TPR(Beta_hat, Beta)

  KL_divergence = mean(sapply(X_list_test, function(X) mean(kl_divergence(compute_probabilities_no_Gamma(X, alpha_hat, Beta_hat),
                                                                          compute_probabilities_no_Gamma(X, alpha, Beta)))))

  hellinger_distance = mean(sapply(X_list_test, function(X) mean(hellinger_distance(compute_probabilities_no_Gamma(X, alpha_hat, Beta_hat),
                                                                                    compute_probabilities_no_Gamma(X, alpha, Beta)))))

  error = mean(sapply(1:length(Y_list_test), function(k) error(predict_categories(X_list_test[[k]], alpha_hat, Beta_hat, categories), Y_list_test[[k]])))

  return(list(Beta_SSE = Beta_SSE, Beta_FPR = Beta_FPR, Beta_TPR = Beta_TPR, KL_divergence = KL_divergence, hellinger_distance = hellinger_distance, error = error, alpha = alpha, Beta = Beta, alpha_hat = alpha_hat, Beta_hat = Beta_hat))

}

compute_best_case_performance = function(alpha, Beta, all_alpha_hats, all_Beta_hats, Y_list_test, X_list_test, categories) {

  keep = c("Beta_SSE", "KL_divergence", "hellinger_distance", "error")

  result = mapply(alpha_hat = all_alpha_hats, Beta_hat = all_Beta_hats, FUN = function(alpha_hat, Beta_hat) {
    compute_performance(alpha, Beta, alpha_hat, Beta_hat, Y_list_test, X_list_test, categories)[keep]
  }, SIMPLIFY = FALSE)

  result = matrix(unlist(result), nrow = length(keep))

  result = apply(result, 1, min)
  names(result) = keep

  return(result)

}

get_all_alphas_Betas_IBMR = function(fit) {

  alphas = list()
  Betas = list()

  for (r in 1:fit$n_rho) {

    for (l in 1:fit$n_lambda) {

      model = fit$model_fits[[r]][[l]]

      if (!is.null(model)) {

        alphas = c(alphas, list(fit$model_fits[[r]][[l]]$alpha))
        Betas = c(Betas, list(fit$model_fits[[r]][[l]]$Beta))

      }

    }

  }

  return(list(all_alphas = alphas, all_Betas = Betas))

}

get_all_alphas_Betas_IBMR_no_Gamma = function(fit) {

  alphas = list()
  Betas = list()

  for (l in 1:fit$n_lambda) {

    model = fit$model_fits[[l]]

    if (!is.null(model)) {

      alphas = c(alphas, list(fit$model_fits[[l]]$alpha))
      Betas = c(Betas, list(fit$model_fits[[l]]$Beta))

    }

  }

  return(list(all_alphas = alphas, all_Betas = Betas))

}

get_all_alphas_Betas_group_lasso = function(fit) {

  alphas = list()
  Betas = list()

  a = 1

  for (l in 1:fit$n_lambda) {

    model = fit$model_fits[[a]][[l]]

    if (!is.null(model)) {

      alphas = c(alphas, list(fit$model_fits[[a]][[l]]$alpha))
      Betas = c(Betas, list(fit$model_fits[[a]][[l]]$Beta))

    }

  }

  return(list(all_alphas = alphas, all_Betas = Betas))

}

prepare_output_IBMR = function(fit) {

  IBMR = list(
    alpha = fit$best_model$alpha,
    Beta = fit$best_model$Beta,
    tuning_parameters = fit$best_tuning_parameters,
    validation_negative_log_likelihood = fit$validation_negative_log_likelihood
  )
  IBMR = c(IBMR, get_all_alphas_Betas_IBMR(fit))

  return(list(IBMR = IBMR))

}

prepare_output_IBMR_no_Gamma = function(fit) {

  IBMR = list(
    alpha = fit$best_model$alpha,
    Beta = fit$best_model$Beta,
    tuning_parameters = fit$best_tuning_parameters,
    validation_negative_log_likelihood = fit$validation_negative_log_likelihood
  )
  IBMR = c(IBMR, get_all_alphas_Betas_IBMR_no_Gamma(fit))

  return(list(IBMR_no_Gamma = IBMR))

}

prepare_output_glmnet = function(fit) {

  elastic_net = list(
    alpha = fit$best_model$alpha,
    Beta = fit$best_model$Beta,
    tuning_parameters = fit$best_tuning_parameters,
    validation_negative_log_likelihood = fit$validation_negative_log_likelihood
  )
  elastic_net = c(elastic_net, get_all_alphas_Betas_IBMR(fit))

  group_lasso = list(
    alpha = fit$best_model_group_lasso$alpha,
    Beta = fit$best_model_group_lasso$Beta,
    tuning_parameters = fit$best_tuning_parameters_group_lasso,
    validation_negative_log_likelihood = fit$validation_negative_log_likelihood[1,]
  )
  group_lasso = c(group_lasso, get_all_alphas_Betas_group_lasso(fit))

  return(list(elastic_net = elastic_net, group_lasso = group_lasso))

}

fit_IBMR = function(data) {

  fit = IBMR(
    data$train$Y_list,
    data$train$category_mappings$categories,
    data$train$category_mappings$category_mappings,
    data$train$X_list,
    data$train$Z_list,
    data$validation$Y_list,
    data$validation$category_mappings$category_mappings,
    data$validation$X_list
  )

  return(prepare_output_IBMR(fit))

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

  return(prepare_output_IBMR_no_Gamma(fit))

}

fit_elastic_net = function(data) {

  fit = fit_glmnet(
    data$train$Y_list,
    data$train$category_mappings$categories,
    data$train$X_list,
    data$validation$Y_list,
    data$validation$X_list
  )

  return(prepare_output_glmnet(fit))

}

fit_ORACLE = function(data) {

  list(
    ORACLE = list(
      alpha = data$alpha,
      Beta = data$Beta,
      tuning_parameters = 0,
      validation_negative_log_likelihood = 0,
      all_alphas = list(data$alpha),
      all_Betas = list(data$Beta)
    )
  )

}
