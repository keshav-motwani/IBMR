generate_simulation_data = function(category_mappings,
                                    N,
                                    p,
                                    nonzero,
                                    b,
                                    rank,
                                    batch_effect,
                                    q = NULL,
                                    K = NULL,
                                    replicate) {

  set.seed(replicate, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")

  if (is.null(category_mappings)) {

    categories = paste0("C", 1:q)
    names(categories) = categories

    category_mappings = create_fine_category_mappings(categories, K)
    category_mappings_fine = category_mappings

  } else {

    K = length(category_mappings$category_mappings)
    category_mappings_fine = create_fine_category_mappings(category_mappings$categories, K)

  }

  alpha = simulate_alpha(category_mappings$categories)
  Beta = simulate_Beta(category_mappings$categories, p, nonzero, -b, b)

  X_star_list = simulate_X_star_list(rep(N / K, K), p)
  Y_list = simulate_Y_list(category_mappings$categories, category_mappings$inverse_category_mappings, X_star_list, alpha, Beta)

  X_star_list_val = simulate_X_star_list(rep(N / K, K), p)
  Y_list_val = simulate_Y_list(category_mappings$categories, category_mappings$inverse_category_mappings, X_star_list_val, alpha, Beta)

  X_star_list_test = simulate_X_star_list(N, p)
  Y_list_test = simulate_Y_list(category_mappings_fine$categories, category_mappings_fine$inverse_category_mappings[1], X_star_list_test, alpha, Beta)

  U_list = simulate_U_list(X_star_list, rank, batch_effect)
  X_list = compute_X_list(X_star_list, U_list)

  U_list_val = simulate_U_list(X_star_list_val, rank, batch_effect)
  X_list_val = compute_X_list(X_star_list_val, U_list_val)

  output = list(
    observed = list(
      train = list(
        X_list = X_list,
        Y_list = Y_list,
        category_mappings = category_mappings
      ),
      validation = list(
        X_list = X_list_val,
        Y_list = Y_list_val,
        category_mappings = category_mappings
      )
    ),
    ORC_clean = list(
      train = list(
        X_list = X_star_list,
        Y_list = Y_list,
        category_mappings = category_mappings
      ),
      validation = list(
        X_list = X_star_list_val,
        Y_list = Y_list_val,
        category_mappings = category_mappings
      )
    ),
    ORC_fine = list(
      train = list(
        X_list = X_list,
        Y_list = get_fine_categories(Y_list),
        category_mappings = category_mappings_fine
      ),
      validation = list(
        X_list = X_list_val,
        Y_list = get_fine_categories(Y_list_val),
        category_mappings = category_mappings_fine
      )
    ),
    ORC_fine_clean = list(
      train = list(
        X_list = X_star_list,
        Y_list = get_fine_categories(Y_list),
        category_mappings = category_mappings_fine
      ),
      validation = list(
        X_list = X_star_list_val,
        Y_list = get_fine_categories(Y_list_val),
        category_mappings = category_mappings_fine
      )
    )
  )

  for (ORC_type in names(output)) {

    output[[ORC_type]]$test = list(X_star_list = X_star_list_test,
                                   Y_list_fine = Y_list_test)

    output[[ORC_type]]$alpha = alpha
    output[[ORC_type]]$Beta = Beta

  }

  return(output)

}

generate_simulation_data_from_real = function(category_mappings,
                                              X_star,
                                              glmnet_fit,
                                              N,
                                              sparsity,
                                              rank,
                                              batch_effect,
                                              replicate) {

  set.seed(replicate, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")

  coef = extract_alpha_Beta_from_glmnet(glmnet_fit, sparsity)
  alpha = coef$alpha
  Beta = coef$Beta

  K = length(category_mappings$category_mappings)

  categories = colnames(Beta)
  category_mappings_fine = create_fine_category_mappings(categories, K)

  X_star = X_star[sample(1:nrow(X_star), nrow(X_star)), ]

  n_k = c(rep(N / K, K), rep(N / K, K), N)
  indices_list = lapply(2:length(n_k), function(i) (sum(n_k[1:i-1]) + 1):sum(n_k[1:i]))
  indices_list = c(list(1:(n_k[1])), indices_list)

  X_star_list = lapply(indices_list[1:K], function(indices) X_star[indices, ])
  Y_list = simulate_Y_list(categories, category_mappings$inverse_category_mappings, X_star_list, alpha, Beta)

  X_star_list_val = lapply(indices_list[(K + 1):(2 * K)], function(indices) X_star[indices, ])
  Y_list_val = simulate_Y_list(categories, category_mappings$inverse_category_mappings, X_star_list_val, alpha, Beta)

  X_star_list_test = lapply(indices_list[2 * K + 1], function(indices) X_star[indices, ])
  Y_list_test = simulate_Y_list(categories, category_mappings$inverse_category_mappings, X_star_list_test, alpha, Beta)

  U_list = simulate_U_list(X_star_list, rank, batch_effect)
  X_list = compute_X_list(X_star_list, U_list)

  U_list_val = simulate_U_list(X_star_list_val, rank, batch_effect)
  X_list_val = compute_X_list(X_star_list_val, U_list_val)

  output = list(
    observed = list(
      train = list(
        X_list = X_list,
        Y_list = Y_list,
        category_mappings = category_mappings
      ),
      validation = list(
        X_list = X_list_val,
        Y_list = Y_list_val,
        category_mappings = category_mappings
      )
    ),
    ORC_clean = list(
      train = list(
        X_list = X_star_list,
        Y_list = Y_list,
        category_mappings = category_mappings
      ),
      validation = list(
        X_list = X_star_list_val,
        Y_list = Y_list_val,
        category_mappings = category_mappings
      )
    ),
    ORC_fine = list(
      train = list(
        X_list = X_list,
        Y_list = get_fine_categories(Y_list),
        category_mappings = category_mappings_fine
      ),
      validation = list(
        X_list = X_list_val,
        Y_list = get_fine_categories(Y_list_val),
        category_mappings = category_mappings_fine
      )
    ),
    ORC_fine_clean = list(
      train = list(
        X_list = X_star_list,
        Y_list = get_fine_categories(Y_list),
        category_mappings = category_mappings_fine
      ),
      validation = list(
        X_list = X_star_list_val,
        Y_list = get_fine_categories(Y_list_val),
        category_mappings = category_mappings_fine
      )
    )
  )

  for (ORC_type in names(output)) {

    output[[ORC_type]]$test = list(X_star_list = X_star_list_test,
                                   Y_list_fine = Y_list_test)

    output[[ORC_type]]$alpha = alpha
    output[[ORC_type]]$Beta = Beta

  }

  return(output)

}

extract_alpha_Beta_from_glmnet = function(glmnet_fit, sparsity) {

  library(glmnet)

  p = glmnet_fit$glmnet.fit$dim[1]

  lambda_index = which.min(abs(glmnet_fit$nzero/p - sparsity))

  coef = as.matrix(do.call(cbind, coef(glmnet_fit, s = glmnet_fit$lambda[lambda_index])))
  colnames(coef) = glmnet_fit$glmnet.fit$classnames

  return(list(alpha = coef[1, , drop = TRUE], Beta = coef[-1, , drop = FALSE]))

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
  ORC_type = paste0("ORC", strsplit(method, "ORC")[[1]][2])
  ORC_type = ifelse(ORC_type == "ORCNA", "observed", ORC_type)
  method = strsplit(strsplit(method, "_ORC")[[1]][1], "fit_")[[1]][1]

  method_function = get(paste0("fit_", method))

  data = do.call(simulation_function, parameters[intersect(names(parameters), formalArgs(simulation_function))])[[ORC_type]]
  data$parameters = parameters

  fits = method_function(data)

  results = vector("list", length(fits))

  for (i in 1:length(fits)) {

    fit = fits[[i]]

    performance = compute_performance(data$test$Y_list_fine,
                                      data$test$X_star_list,
                                      data$alpha,
                                      data$Beta,
                                      fit$alpha_hat,
                                      fit$Beta_hat,
                                      fit$test_estimated_probabilities)

    best_case_performance = compute_best_case_performance(data$test$Y_list_fine,
                                                          data$test$X_star_list,
                                                          data$alpha,
                                                          data$Beta,
                                                          fit$all_alpha_hats,
                                                          fit$all_Beta_hats,
                                                          fit$all_test_estimated_probabilities)

    fit$test_estimated_probabilities = NULL

    fit$all_alpha_hats = NULL
    fit$all_Beta_hats = NULL
    fit$all_test_estimated_probabilities = NULL

    if (length(fits) > 1) parameters$method = paste0(parameters$method, "_", names(fits)[i])

    results[[i]] = list(parameters = parameters, performance = performance, best_case_performance = best_case_performance, fit = fit, true = list(alpha = data$alpha, Beta = data$Beta))

  }

  return(results)

}

compute_performance = function(Y_list_test, X_list_test, alpha, Beta, alpha_hat, Beta_hat, test_estimated_probabilities) {

  if (!is.null(alpha_hat)) {

    Beta_SSE = Beta_SSE(Beta_hat, Beta)
    Beta_FPR = Beta_FPR(Beta_hat, Beta)
    Beta_TPR = Beta_TPR(Beta_hat, Beta)

  } else {

    Beta_SSE = Beta_FPR = Beta_TPR = NA

  }

  KL_divergence = mean(unlist(mapply(test_estimated_probabilities, X_list_test, FUN = function(P_hat, X) mean(kl_divergence(P_hat, IBMR:::compute_probabilities_no_Gamma(X, alpha, Beta))), SIMPLIFY = FALSE)))

  hellinger_distance = mean(unlist(mapply(test_estimated_probabilities, X_list_test, FUN = function(P_hat, X) mean(hellinger_distance(P_hat, IBMR:::compute_probabilities_no_Gamma(X, alpha, Beta))), SIMPLIFY = FALSE)))

  predicted_categories = predict_categories(test_estimated_probabilities)

  error = mean(unlist(mapply(predicted_categories, Y_list_test, FUN = function(predictions, Y) error(predictions, Y), SIMPLIFY = FALSE)))

  return(list(Beta_SSE = Beta_SSE, Beta_FPR = Beta_FPR, Beta_TPR = Beta_TPR, KL_divergence = KL_divergence, hellinger_distance = hellinger_distance, error = error))

}

compute_best_case_performance = function(Y_list_test, X_list_test, alpha, Beta, all_alpha_hats, all_Beta_hats, all_test_estimated_probabilities) {

  keep = c("Beta_SSE", "KL_divergence", "hellinger_distance", "error")

  if (is.null(all_alpha_hats) & is.null(all_Beta_hats)) all_alpha_hats = all_Beta_hats = vector("list", length(all_test_estimated_probabilities))

  result = mapply(alpha_hat = all_alpha_hats, Beta_hat = all_Beta_hats, test_estimated_probabilities = all_test_estimated_probabilities, FUN = function(alpha_hat, Beta_hat, test_estimated_probabilities) {
    compute_performance(Y_list_test, X_list_test, alpha, Beta, alpha_hat, Beta_hat, test_estimated_probabilities)[keep]
  }, SIMPLIFY = FALSE)

  result = matrix(unlist(result), nrow = length(keep))

  result = apply(result, 1, min)
  names(result) = keep

  return(result)

}

get_all_estimates_two_tuning_parameters = function(fit, X_list_test) {

  alphas = list()
  Betas = list()
  probabilities = list()

  for (r in 1:fit$n_rho) {

    for (l in 1:fit$n_lambda) {

      model = fit$model_fits[[r]][[l]]

      if (!is.null(model)) {

        alphas = c(alphas, list(model$alpha))
        Betas = c(Betas, list(model$Beta))
        probabilities = c(probabilities, list(predict_probabilities(model, X_list_test)))

      }

    }

  }

  return(list(all_alpha_hats = alphas, all_Beta_hats = Betas, all_test_estimated_probabilities = probabilities))

}

get_all_estimates_one_tuning_parameter = function(fit, X_list_test) {

  alphas = list()
  Betas = list()
  probabilities = list()

  for (l in 1:fit$n_lambda) {

    model = fit$model_fits[[l]]

    if (!is.null(model)) {

      alphas = c(alphas, list(model$alpha))
      Betas = c(Betas, list(model$Beta))
      probabilities = c(probabilities, list(predict_probabilities(model, X_list_test)))

    }

  }

  return(list(all_alpha_hats = alphas, all_Beta_hats = Betas, all_test_estimated_probabilities = probabilities))

}

prepare_output_IBMR = function(fit, X_list_test) {

  IBMR = list(
    alpha_hat = fit$best_model$alpha,
    Beta_hat = fit$best_model$Beta,
    test_estimated_probabilities = predict_probabilities(fit$best_model, X_list_test),
    tuning_parameters = fit$best_tuning_parameters,
    validation_negative_log_likelihood = fit$validation_negative_log_likelihood,
    best_model = fit$best_model
  )
  IBMR = c(IBMR, get_all_estimates_two_tuning_parameters(fit, X_list_test))

  return(list(IBMR = IBMR))

}

prepare_output_IBMR_no_Gamma = function(fit, X_list_test) {

  IBMR = list(
    alpha_hat = fit$best_model$alpha,
    Beta_hat = fit$best_model$Beta,
    test_estimated_probabilities = predict_probabilities(fit$best_model, X_list_test),
    tuning_parameters = fit$best_tuning_parameters,
    validation_negative_log_likelihood = fit$validation_negative_log_likelihood,
    best_model = fit$best_model
  )
  IBMR = c(IBMR, get_all_estimates_one_tuning_parameter(fit, X_list_test))

  return(list(IBMR_no_Gamma = IBMR))

}

prepare_output_glmnet_split = function(fit, X_list_test) {

  glmnet_split = list(
    test_estimated_probabilities = predict_probabilities_glmnet_split(fit, X_list_test),
    all_test_estimated_probabilities = list(predict_probabilities_glmnet_split(fit, X_list_test)),
    best_model = fit
  )

  return(list(glmnet_split = glmnet_split))

}

fit_IBMR = function(data) {

  fit = IBMR(
    data$train$Y_list,
    data$train$category_mappings$categories,
    data$train$category_mappings$category_mappings,
    data$train$X_list,
    compute_pca_for_Z_list(data$train$X_list, ifelse(data$parameters$p > 100, 50, 10), TRUE),
    data$validation$Y_list,
    data$validation$category_mappings$category_mappings,
    data$validation$X_list
  )

  return(prepare_output_IBMR(fit, data$test$X_star_list))

}

fit_IBMR_int = function(data) {

  fit = IBMR(
    data$train$Y_list,
    data$train$category_mappings$categories,
    data$train$category_mappings$category_mappings,
    data$train$X_list,
    lapply(data$train$X_list, function(X) matrix(1, nrow = nrow(X), ncol = 1)),
    data$validation$Y_list,
    data$validation$category_mappings$category_mappings,
    data$validation$X_list
  )

  return(prepare_output_IBMR(fit, data$test$X_star_list))

}

fit_IBMR_common_Gamma = function(data) {

  fit = IBMR(
    data$train$Y_list,
    data$train$category_mappings$categories,
    data$train$category_mappings$category_mappings,
    data$train$X_list,
    compute_pca_for_Z_list(data$train$X_list, ifelse(data$parameters$p > 100, 50, 10), TRUE),
    data$validation$Y_list,
    data$validation$category_mappings$category_mappings,
    data$validation$X_list,
    common_Gamma = TRUE
  )

  return(prepare_output_IBMR(fit, data$test$X_star_list))

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

  return(prepare_output_IBMR_no_Gamma(fit, data$test$X_star_list))

}

fit_glmnet_subset = function(data) {

  fit = glmnet_subset(
    data$train$Y_list,
    data$train$category_mappings$categories,
    data$train$category_mappings$category_mappings,
    data$train$X_list,
    data$validation$Y_list,
    data$validation$category_mappings$category_mappings,
    data$validation$X_list
  )

  return(prepare_output_IBMR_no_Gamma(fit, data$test$X_star_list))

}

fit_glmnet_split = function(data) {

  fit = glmnet_split(
    data$train$Y_list,
    data$train$category_mappings$categories,
    data$train$category_mappings$category_mappings,
    data$train$X_list,
    data$validation$Y_list,
    data$validation$category_mappings$category_mappings,
    data$validation$X_list
  )

  return(prepare_output_glmnet_split(fit, data$test$X_star_list))

}

fit_glmnet_relabel = function(data) {

  fit = glmnet_relabel(
    data$train$Y_list,
    data$train$category_mappings$categories,
    data$train$category_mappings$category_mappings,
    data$train$X_list,
    data$validation$Y_list,
    data$validation$category_mappings$category_mappings,
    data$validation$X_list
  )

  return(prepare_output_IBMR(fit, data$test$X_star_list))

}

fit_ORACLE = function(data) {

  list(
    ORACLE = list(
      alpha_hat = data$alpha,
      Beta_hat = data$Beta,
      test_estimated_probabilities = predict_probabilities(list(alpha = data$alpha, Beta = data$Beta), data$test$X_star_list),
      tuning_parameters = 0,
      validation_negative_log_likelihood = 0,
      all_alphas = list(data$alpha),
      all_Betas = list(data$Beta),
      all_test_estimated_probabilities = list(predict_probabilities(list(alpha = data$alpha, Beta = data$Beta), data$test$X_star_list))
    )
  )

}
