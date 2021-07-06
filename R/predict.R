#' @export
predict_probabilities = function(model, X_list) {

  alpha = model$alpha
  Beta = model$Beta

  categories = colnames(Beta)

  probabilities = vector("list", length(X_list))

  for (k in 1:length(predictions)) {

    P = compute_probabilities_no_Gamma(X_list[[k]], alpha, Beta)
    colnames(P) = categories

    probabilities[[k]] = P

  }

  return(probabilities)

}

#' @export
predict_conditional_probabilities = function(model, Y_list, category_mappings, X_list) {

  alpha = model$alpha
  Beta = model$Beta

  categories = colnames(Beta)

  conditional_probabilities = vector("list", length(X_list))

  for (k in 1:length(predictions)) {

    P = compute_probabilities_no_Gamma(X_list[[k]], alpha, Beta)
    Y_matrix = create_Y_matrix(Y_list[[k]], categories, category_mappings[[k]])
    C = compute_conditional_probabilities(Y_matrix, P)
    colnames(C) = categories

    conditional_probabilities[[k]] = C

  }

  return(conditional_probabilities)

}

#' @export
predict_probabilities_train = function(model, X_list, Z_list) {

  alpha = model$alpha
  Beta = model$Beta
  Gamma_list = model$Gamma_list

  categories = colnames(Beta)

  Z_list = standardize_Z(Z_list)

  probabilities = vector("list", length(X_list))

  for (k in 1:length(predictions)) {

    P = compute_probabilities(X_list[[k]], Z_list[[k]], alpha, Beta, Gamma_list[[k]])
    colnames(P) = categories

    probabilities[[k]] = P

  }

  return(probabilities)

}

#' @export
predict_conditional_probabilities_train = function(model, Y_list, category_mappings, X_list, Z_list) {

  alpha = model$alpha
  Beta = model$Beta
  Gamma_list = model$Gamma_list

  stopifnot(length(Gamma_list) == length(X_list))

  categories = colnames(Beta)

  Z_list = standardize_Z(Z_list)

  conditional_probabilities = vector("list", length(X_list))

  for (k in 1:length(predictions)) {

    P = compute_probabilities(X_list[[k]], Z_list[[k]], alpha, Beta, Gamma_list[[k]])
    Y_matrix = create_Y_matrix(Y_list[[k]], categories, category_mappings[[k]])
    C = compute_conditional_probabilities(Y_matrix, P)
    colnames(C) = categories

    conditional_probabilities[[k]] = C

  }

  return(conditional_probabilities)

}

#' @export
predict_categories = function(predicted_probabilities) {

    categories = colnames(predicted_probabilities[[1]])

    predictions = vector("list", length(predicted_probabilities))

    for (k in 1:length(predictions)) {

      predictions[[k]] = categories[apply(predicted_probabilities[[k]], 1, which.max)]

    }

    return(predictions)

}
