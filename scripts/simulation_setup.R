generate_data_random_X_and_structured_Beta = function(category_mappings,
                                                      N,
                                                      p,
                                                      nonzero,
                                                      s,
                                                      sigma,
                                                      rank,
                                                      batch_effect,
                                                      replicate) {

  set.seed(replicate, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")

  alpha = simulate_alpha(category_mappings$categories)
  Beta = simulate_structured_Beta(category_mappings$splits_per_level, p, nonzero, c(s, nonzero - s), sigma)

  K = length(category_mappings$category_mappings)
  category_mappings_fine = create_fine_category_mappings(category_mappings$categories, K)

  X_star_list = simulate_X_star_list(rep(N / K, K), p)
  Y_list = simulate_Y_list(category_mappings$categories, category_mappings$inverse_category_mappings, X_star_list, alpha, Beta)

  X_star_list_val = simulate_X_star_list(rep(N / K, K), p)
  Y_list_val = simulate_Y_list(category_mappings$categories, category_mappings$inverse_category_mappings, X_star_list_val, alpha, Beta)

  X_star_list_test = simulate_X_star_list(10000, p)
  Y_list_test = simulate_Y_list(category_mappings_fine$categories, category_mappings_fine$inverse_category_mappings[1], X_star_list_test, alpha, Beta)

  U_list = simulate_U_list(X_star_list, rank, batch_effect)
  X_list = compute_X_list(X_star_list, U_list)

  U_list_val = simulate_U_list(X_star_list_val, rank, batch_effect)
  X_list_val = compute_X_list(X_star_list_val, U_list_val)

  output = prepare_data(Y_list = Y_list,
                        category_mappings = category_mappings,
                        category_mappings_fine = category_mappings_fine,
                        X_list = X_list,
                        X_star_list = X_star_list,
                        Y_list_validation = Y_list_val,
                        category_mappings_validation = category_mappings,
                        category_mappings_fine_validation = category_mappings_fine,
                        X_list_validation = X_list_val,
                        X_star_list_validation = X_star_list_val,
                        Y_list_test = Y_list_test,
                        category_mappings_test = list(categories = category_mappings_fine$categories, category_mappings = category_mappings_fine$category_mappings[1], inverse_category_mappings = category_mappings_fine$inverse_category_mappings[1]),
                        X_list_test = X_star_list_test,
                        alpha = alpha,
                        Beta = Beta)

  return(output)

}

prepare_data = function(Y_list,
                        category_mappings,
                        category_mappings_fine,
                        X_list,
                        X_star_list,
                        sample_list = NULL,
                        Y_list_validation,
                        category_mappings_validation,
                        category_mappings_fine_validation,
                        X_list_validation,
                        X_star_list_validation,
                        Y_list_test,
                        category_mappings_test,
                        X_list_test,
                        alpha = NULL,
                        Beta = NULL) {

  sapply(Y_list, function(x) print(table(x)))
  print(table(unlist((Y_list))))
  print(table(unlist(get_fine_categories(Y_list))))

  output = list(
    observed = list(
      train = list(
        X_list = X_list,
        Y_list = Y_list,
        category_mappings = category_mappings,
        sample_list = sample_list
      ),
      validation = list(
        X_list = X_list_validation,
        Y_list = Y_list_validation,
        category_mappings = category_mappings_validation
      )
    ),
    ORC_clean = list(
      train = list(
        X_list = X_star_list,
        Y_list = Y_list,
        category_mappings = category_mappings
      ),
      validation = list(
        X_list = X_star_list_validation,
        Y_list = Y_list_validation,
        category_mappings = category_mappings_validation
      )
    ),
    ORC_fine = list(
      train = list(
        X_list = X_list,
        Y_list = get_fine_categories(Y_list),
        category_mappings = category_mappings_fine
      ),
      validation = list(
        X_list = X_list_validation,
        Y_list = get_fine_categories(Y_list_validation),
        category_mappings = category_mappings_fine_validation
      )
    ),
    ORC_fine_clean = list(
      train = list(
        X_list = X_star_list,
        Y_list = get_fine_categories(Y_list),
        category_mappings = category_mappings_fine
      ),
      validation = list(
        X_list = X_star_list_validation,
        Y_list = get_fine_categories(Y_list_validation),
        category_mappings = category_mappings_fine_validation
      )
    )
  )

  for (ORC_type in names(output)) {

    output[[ORC_type]]$test = list(X_list = X_list_test,
                                   Y_list = Y_list_test,
                                   category_mappings = category_mappings_test)

    output[[ORC_type]]$alpha = alpha
    output[[ORC_type]]$Beta = Beta

  }

  return(output)

}

extract_alpha_Beta_from_glmnet = function(glmnet_fit, nonsparsity) {

  library(glmnet)

  p = nrow(glmnet_fit$beta[[1]])

  lambda_index = which.min(abs(glmnet_fit$df/p - nonsparsity))[1]

  stopifnot(abs(glmnet_fit$df[lambda_index]/p - nonsparsity) < 0.05)

  coef = as.matrix(do.call(cbind, coef(glmnet_fit, s = glmnet_fit$lambda[lambda_index])))
  colnames(coef) = glmnet_fit$classnames

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

  print("Preparing data")

  data = do.call(simulation_function, parameters[intersect(names(parameters), formalArgs(simulation_function))])[[ORC_type]]
  data$parameters = parameters

  print("Fitting model")

  time = system.time({fits = method_function(data)})[3]

  results = vector("list", length(fits))

  print("Evaluating performance")

  for (i in 1:length(fits)) {

    fit = fits[[i]]

    performance = compute_performance(data$test$Y_list,
                                      data$test$category_mappings,
                                      data$test$X_list,
                                      data$alpha,
                                      data$Beta,
                                      fit$alpha_hat,
                                      fit$Beta_hat,
                                      fit$test_estimated_probabilities)
    performance$time = time

    # best_case_performance = compute_best_case_performance(data$test$Y_list,
    #                                                       data$test$category_mappings,
    #                                                       data$test$X_list,
    #                                                       data$alpha,
    #                                                       data$Beta,
    #                                                       fit$all_alpha_hats,
    #                                                       fit$all_Beta_hats,
    #                                                       fit$all_test_estimated_probabilities)

    fit$test_estimated_probabilities = NULL

    fit$all_alpha_hats = NULL
    fit$all_Beta_hats = NULL
    fit$all_test_estimated_probabilities = NULL

    if (length(fits) > 1) parameters$method = paste0(parameters$method, "_", names(fits)[i])

    parameters$X_star = NULL
    parameters$glmnet_fit = NULL

    results[[i]] = list(parameters = parameters, performance = performance, best_case_performance = NULL, fit = fit, true = list(alpha = data$alpha, Beta = data$Beta))

  }

  return(results)

}

compute_performance = function(Y_list_test, category_mappings_test, X_list_test, alpha, Beta, alpha_hat, Beta_hat, test_estimated_probabilities) {

  if (!is.null(Beta)) {

    if (!is.null(Beta_hat)) {

      Beta_SSE = Beta_SSE(Beta_hat, Beta)
      Beta_FPR = Beta_FPR(Beta_hat, Beta)
      Beta_TPR = Beta_TPR(Beta_hat, Beta)

    } else {

      Beta_SSE = Beta_FPR = Beta_TPR = NA

    }

    KL_divergence = mean(unlist(mapply(test_estimated_probabilities, X_list_test, FUN = function(P_hat, X) mean(kl_divergence(P_hat, IBMR:::compute_probabilities_no_Gamma(X, alpha, Beta))), SIMPLIFY = FALSE)))

    hellinger_distance = mean(unlist(mapply(test_estimated_probabilities, X_list_test, FUN = function(P_hat, X) mean(hellinger_distance(P_hat, IBMR:::compute_probabilities_no_Gamma(X, alpha, Beta))), SIMPLIFY = FALSE)))

  } else {

    Beta_SSE = Beta_FPR = Beta_TPR = KL_divergence = hellinger_distance = NA

  }

  predicted_categories = predict_categories(test_estimated_probabilities, category_mappings_test$category_mappings)

  error = error(unlist(predicted_categories), unlist(Y_list_test))
  balanced_error = balanced_error(unlist(predicted_categories), unlist(Y_list_test))

  Y_matrix_list_test = mapply(Y = Y_list_test, category_mapping = category_mappings_test$category_mappings, FUN = function(Y, category_mapping) create_Y_matrix(Y, category_mappings_test$categories, category_mapping), SIMPLIFY = FALSE)
  # nll = compute_negative_log_likelihood_no_Gamma(Y_matrix_list_test, X_list_test, alpha_hat, Beta_hat, length(unlist(Y_list_test)))
  nll = compute_negative_log_likelihood_from_probabilities(Y_matrix_list_test, test_estimated_probabilities, length(unlist(Y_list_test)))

  confusion_matrix = table(unlist(Y_list_test), unlist(predicted_categories))

  return(list(Beta_SSE = Beta_SSE, Beta_FPR = Beta_FPR, Beta_TPR = Beta_TPR, KL_divergence = KL_divergence, hellinger_distance = hellinger_distance, error = error, balanced_error = balanced_error, nll = nll, confusion_matrix = confusion_matrix))

}

compute_best_case_performance = function(Y_list_test, category_mappings_test, X_list_test, alpha, Beta, all_alpha_hats, all_Beta_hats, all_test_estimated_probabilities) {

  keep = c("Beta_SSE", "KL_divergence", "hellinger_distance", "error", "balanced_error", "nll")

  if (is.null(all_alpha_hats) & is.null(all_Beta_hats)) all_alpha_hats = all_Beta_hats = vector("list", length(all_test_estimated_probabilities))

  result = mapply(alpha_hat = all_alpha_hats, Beta_hat = all_Beta_hats, test_estimated_probabilities = all_test_estimated_probabilities, FUN = function(alpha_hat, Beta_hat, test_estimated_probabilities) {
    compute_performance(Y_list_test, category_mappings_test, X_list_test, alpha, Beta, alpha_hat, Beta_hat, test_estimated_probabilities)[keep]
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
    X_mean = fit$X_mean, X_sd = fit$X_sd,
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
    X_mean = fit$X_mean, X_sd = fit$X_sd,
    best_model = fit$best_model
  )
  IBMR = c(IBMR, get_all_estimates_one_tuning_parameter(fit, X_list_test))

  return(list(IBMR_no_Gamma = IBMR))

}

prepare_output_relabel = function(fit, X_list_test) {

  return(list(subset = prepare_output_IBMR_no_Gamma(attr(fit, "subset")), relabel = prepare_output_IBMR_no_Gamma(fit)))

}

fit_IBMR = function(data) {

  fit = IBMR(
    data$train$Y_list,
    data$train$category_mappings$categories,
    data$train$category_mappings$category_mappings,
    data$train$X_list,
    compute_pca_for_Z_list(data$train$X_list, ifelse(ncol(data$train$X_list[[1]]) > 100, 50, 10), TRUE),
    data$validation$Y_list,
    data$validation$category_mappings$category_mappings,
    data$validation$X_list
  )

  return(prepare_output_IBMR(fit, data$test$X_list))

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

  return(prepare_output_IBMR(fit, data$test$X_list))

}

fit_IBMR_int_1en6 = function(data) {

  fit = IBMR(
    data$train$Y_list,
    data$train$category_mappings$categories,
    data$train$category_mappings$category_mappings,
    data$train$X_list,
    lapply(data$train$X_list, function(X) matrix(1, nrow = nrow(X), ncol = 1)),
    data$validation$Y_list,
    data$validation$category_mappings$category_mappings,
    data$validation$X_list,
    rho_min_ratio = 1e-6
  )

  return(prepare_output_IBMR(fit, data$test$X_list))

}

fit_IBMR_common_Gamma = function(data) {

  fit = IBMR(
    data$train$Y_list,
    data$train$category_mappings$categories,
    data$train$category_mappings$category_mappings,
    data$train$X_list,
    compute_pca_for_Z_list(data$train$X_list, ifelse(ncol(data$train$X_list[[1]]) > 100, 50, 10), TRUE),
    data$validation$Y_list,
    data$validation$category_mappings$category_mappings,
    data$validation$X_list,
    common_Gamma = TRUE
  )

  return(prepare_output_IBMR(fit, data$test$X_list))

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

  return(prepare_output_IBMR_no_Gamma(fit, data$test$X_list))

}

fit_subset = function(data) {

  fit = IBMR_no_Gamma_subset(
    data$train$Y_list,
    data$train$category_mappings$categories,
    data$train$category_mappings$category_mappings,
    data$train$X_list,
    data$validation$Y_list,
    data$validation$category_mappings$category_mappings,
    data$validation$X_list
  )

  return(prepare_output_IBMR_no_Gamma(fit, data$test$X_list))

}

fit_relabel = function(data) {

  fit = IBMR_no_Gamma_relabel(
    data$train$Y_list,
    data$train$category_mappings$categories,
    data$train$category_mappings$category_mappings,
    data$train$X_list,
    data$validation$Y_list,
    data$validation$category_mappings$category_mappings,
    data$validation$X_list
  )

  return(prepare_output_IBMR_no_Gamma(fit, data$test$X_list))

}

fit_Seurat = function(data) {

  library(Seurat)

  categories = data$train$category_mappings$categories

  subsetted = subset_helper(data$train$Y_list,
                            categories,
                            data$train$category_mappings$category_mappings,
                            data$train$X_list,
                            combine = FALSE)

  data$train$Y_list = subsetted$Y_list
  data$train$X_list = subsetted$X_list

  prepare_Seurat = function(data, which) {
    lapply(1:length(data[[which]]$Y_list), function(k) {
      rownames(data[[which]]$X_list[[k]]) = paste0(sapply(1:nrow(data[[which]]$X_list[[k]]), function(i) paste0(sample(LETTERS, 10), collapse = "")), rownames(data[[which]]$X_list[[k]]))
      seurat_object = CreateSeuratObject(t(data[[which]]$X_list[[k]]))
      seurat_object@assays$RNA@data = t(data[[which]]$X_list[[k]])
      seurat_object$cell_type = data[[which]]$Y_list[[k]]
      seurat_object
    })
  }

  train_datasets = prepare_Seurat(data, "train")
  validation_datasets = prepare_Seurat(data, "validation")
  test_datasets = prepare_Seurat(data, "test")

  features = rownames(train_datasets[[1]])

  n.dim_sequence = 1:5 * 10
  k.anchor_sequence = c(3, 5, 10, 15, 20)

  validation_error = matrix(nrow = length(n.dim_sequence), ncol = length(k.anchor_sequence))

  for (d in 1:length(n.dim_sequence)) {

    for (a in 1:length(k.anchor_sequence)) {

      n.dim = n.dim_sequence[d]
      k.anchor = k.anchor_sequence[a]

      anchors = FindIntegrationAnchors(object.list = train_datasets, dims = 1:n.dim, k.anchor = k.anchor, anchor.features = features)
      integrated = IntegrateData(anchorset = anchors, dims = 1:n.dim)

      DefaultAssay(integrated) = "integrated"
      integrated = ScaleData(integrated, verbose = FALSE, features = features)
      integrated = RunPCA(integrated, npcs = n.dim, verbose = FALSE, features = features)

      P_list_validation = list()

      for (k in 1:length(validation_datasets)) {

        anchors = FindTransferAnchors(reference = integrated, query = validation_datasets[[k]], k.anchor = k.anchor,
                                      dims = 1:n.dim, reference.reduction = "pca", features = features)
        predictions = TransferData(anchorset = anchors, refdata = integrated$cell_type,
                                   dims = 1:n.dim)$predicted.id
        P = create_Y_matrix(predictions, categories, as.list(setNames(categories, categories)))
        P_list_validation = c(P_list_validation, list(P))

      }

      predicted_categories = predict_categories(P_list_validation, data$validation$category_mappings$category_mappings)
      validation_error[d, a] = error(unlist(predicted_categories), unlist(data$validation$Y_list))

      if (validation_error[d, a] <= min(validation_error, na.rm = TRUE)) {
        integrated_best = integrated
        n.dim_best = n.dim
        k.anchor_best = k.anchor
      }

    }

  }

  P_list_test = list()

  for (k in 1:length(test_datasets)) {

    anchors = FindTransferAnchors(reference = integrated_best, query = test_datasets[[k]], k.anchor = k.anchor_best,
                                  dims = 1:n.dim_best, reference.reduction = "pca", features = features)
    predictions = TransferData(anchorset = anchors, refdata = integrated_best$cell_type,
                               dims = 1:n.dim_best)$predicted.id
    P = create_Y_matrix(predictions, categories, as.list(setNames(categories, categories)))
    P_list_test = c(P_list_test, list(P))

  }

  Seurat = list(
    alpha_hat = NULL,
    Beta_hat = NULL,
    test_estimated_probabilities = P_list_test,
    tuning_parameters = which_min(validation_error)[1, ],
    validation_error = validation_error,
    X_mean = NULL, X_sd = NULL,
    best_model = NULL
  )

  return(list(Seurat = Seurat))

}

fit_SingleR = function(data) {

  library(SingleR)

  categories = data$train$category_mappings$categories

  subsetted = subset_helper(data$train$Y_list,
                            categories,
                            data$train$category_mappings$category_mappings,
                            data$train$X_list,
                            combine = FALSE)

  de.n_sequence = c(20, 40, 60, 80, 100)
  quantile_sequence = c(0.6, 0.7, 0.8, 0.9, 1)

  validation_error = matrix(nrow = length(de.n_sequence), ncol = length(quantile_sequence))

  for (d in 1:length(de.n_sequence)) {

    for (q in 1:length(quantile_sequence)) {

      de.n = de.n_sequence[d]
      quantile = quantile_sequence[q]

      fit = trainSingleR(ref = lapply(subsetted$X_list, t), labels = subsetted$Y_list, de.method = "wilcox", aggr.ref = TRUE, de.n = de.n)

      P_list_validation = list()

      for (k in 1:length(data$validation$Y_list)) {

        predictions = classifySingleR(t(data$validation$X_list[[k]]), fit, quantile = quantile)$labels
        P = create_Y_matrix(predictions, categories, as.list(setNames(categories, categories)))
        P_list_validation = c(P_list_validation, list(P))

      }

      predicted_categories = predict_categories(P_list_validation, data$validation$category_mappings$category_mappings)
      validation_error[d, q] = error(unlist(predicted_categories), unlist(data$validation$Y_list))

      if (validation_error[d, q] <= min(validation_error, na.rm = TRUE)) {
        fit_best = fit
        quantile_best = quantile
      }

    }

  }

  P_list_test = list()

  for (k in 1:length(data$test$Y_list)) {

    predictions = classifySingleR(t(data$test$X_list[[k]]), fit_best, quantile = quantile_best)$labels
    P = create_Y_matrix(predictions, categories, as.list(setNames(categories, categories)))
    P_list_test = c(P_list_test, list(P))

  }

  SingleR = list(
    alpha_hat = NULL,
    Beta_hat = NULL,
    test_estimated_probabilities = P_list_test,
    tuning_parameters = which_min(validation_error)[1, ],
    validation_error = validation_error,
    X_mean = NULL, X_sd = NULL,
    best_model = NULL
  )

  return(list(SingleR = SingleR))

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

  return(prepare_output_IBMR_no_Gamma(fit, data$test$X_list))

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

  return(prepare_output_IBMR_no_Gamma(fit, data$test$X_list))

}

fit_ORACLE = function(data) {

  list(
    ORACLE = list(
      alpha_hat = data$alpha,
      Beta_hat = data$Beta,
      test_estimated_probabilities = predict_probabilities(list(alpha = data$alpha, Beta = data$Beta), data$test$X_list),
      tuning_parameters = 0,
      validation_negative_log_likelihood = 0,
      all_alphas = list(data$alpha),
      all_Betas = list(data$Beta),
      all_test_estimated_probabilities = list(predict_probabilities(list(alpha = data$alpha, Beta = data$Beta), data$test$X_list))
    )
  )

}
