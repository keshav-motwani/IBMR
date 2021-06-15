standardize_X = function(X_list) {

  X = do.call(rbind, X_list)
  k = rep(1:length(X_list), sapply(X_list, nrow))

  means = colMeans(X)
  vars = matrixStats::colVars(as.matrix(X)) * (nrow(X) - 1) / nrow(X)

  X = X - tcrossprod(rep(1, nrow(X)), means)
  X = X %*% diag(1 / sqrt(vars))
  colnames(X) = colnames(X_list[[1]])

  X_list = lapply(1:length(X_list), function(i) X[k == i, ])

  attributes(X_list)$mean = means
  attributes(X_list)$sd = sqrt(vars)

  return(X_list)

}

standardize_Z = function(Z_list) {

  means_list = vector("list", length(Z_list))
  sd_list = vector("list", length(Z_list))

  for (k in 1:length(Z_list)) {

    Z = Z_list[[k]]

    means = colMeans(Z)
    vars = matrixStats::colVars(as.matrix(Z)) * (nrow(Z) - 1) / nrow(Z)

    means[vars == 0] = 0
    vars[vars == 0] = 1

    Z = Z - tcrossprod(rep(1, nrow(Z)), means)
    Z = Z %*% diag(1 / sqrt(vars))
    colnames(Z) = colnames(Z_list[[k]])

    means_list[[k]] = means
    sd_list[[k]] = sqrt(vars)
    Z_list[[k]] = Z

  }

  attributes(Z_list)$mean = means_list
  attributes(Z_list)$sd = sd_list

  return(Z_list)

}

adjust_alpha = function(alpha, Beta, X_mean, X_sd) {

  c(alpha - crossprod(Beta, X_mean / X_sd))

}

adjust_Beta = function(Beta, X_sd) {

  Beta = diag(1 / X_sd) %*% Beta

  return(Beta)

}

adjust_fit = function(fit, categories, features, X_mean, X_sd) {

  for (r in 1:fit$n_rho) {
    for (l in 1:fit$n_lambda) {
      if (!is.null(fit$model_fits[[r]][[l]])) {
        fit$model_fits[[r]][[l]]$alpha = adjust_alpha(fit$model_fits[[r]][[l]]$alpha, fit$model_fits[[r]][[l]]$Beta, X_mean, X_sd)
        fit$model_fits[[r]][[l]]$Beta = adjust_Beta(fit$model_fits[[r]][[l]]$Beta, X_sd)
        names(fit$model_fits[[r]][[l]]$alpha) = categories
        colnames(fit$model_fits[[r]][[l]]$Beta) = categories
        rownames(fit$model_fits[[r]][[l]]$Beta) = features
      }
    }
  }

  return(fit)

}

adjust_fit_no_Gamma = function(fit, categories, features, X_mean, X_sd) {

  for (l in 1:fit$n_lambda) {
    fit$model_fits[[l]]$alpha = adjust_alpha(fit$model_fits[[l]]$alpha, fit$model_fits[[l]]$Beta, X_mean, X_sd)
    fit$model_fits[[l]]$Beta = adjust_Beta(fit$model_fits[[l]]$Beta, X_sd)
    names(fit$model_fits[[l]]$alpha) = categories
    colnames(fit$model_fits[[l]]$Beta) = categories
    rownames(fit$model_fits[[l]]$Beta) = features
  }

  return(fit)

}
