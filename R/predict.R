#' @export
predict_probabilities = function(model, X_list) {

  alpha = model$alpha
  Beta = model$Beta

  categories = colnames(Beta)

  probabilities = vector("list", length(X_list))

  for (k in 1:length(probabilities)) {

    P = compute_probabilities_no_Gamma(X_list[[k]], alpha, Beta)
    colnames(P) = categories

    probabilities[[k]] = P

  }

  return(probabilities)

}

#' @export
predict_probabilities_train = function(model, X_list, Z_list) {

  alpha = model$alpha
  Beta = model$Beta
  Gamma_list = model$Gamma_list

  if (length(Gamma_list) == 1) {
    Gamma_list = replicate(length(Y_list), Gamma_list[[1]], simplify = FALSE)
    warning("Assuming that common_Gamma model was fit, if not, make sure that X_list and Z_list that were used for training are used here.")
  }

  if (length(Gamma_list) != length(X_list)) {
    stop("Training data was not provided; make sure that X_list and Z_list that were used for training are used here.")
  }

  stopifnot(length(Gamma_list) == length(X_list))

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
predict_conditional_probabilities = function(predicted_probabilities, Y_list, category_mappings) {

  conditional_probabilities = vector("list", length(Y_list))

  for (k in 1:length(conditional_probabilities)) {

    P = predicted_probabilities[[k]]
    Y_matrix = create_Y_matrix(Y_list[[k]], colnames(P), category_mappings[[k]])
    C = compute_conditional_probabilities(Y_matrix, P)
    colnames(C) = colnames(P)

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
