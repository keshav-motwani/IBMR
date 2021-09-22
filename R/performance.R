#' @export
error = function(predicted, true) {

  mean(true != predicted)

}

#' @export
balanced_error = function(predicted, true) {

  categories = sort(unique(true))
  error = numeric(length(categories))

  for (i in 1:length(categories)) {

    indices = which(true == categories[i])
    true_subsetted = true[indices]
    predicted_subsetted = predicted[indices]

    error[i] = mean(true_subsetted != predicted_subsetted)

  }

  mean(error)

}

#' @export
kl_divergence = function(estimated, true) {

  rowSums(estimated * log(estimated/true))

}

#' @export
hellinger_distance = function(estimated, true) {

  apply(sqrt(estimated) - sqrt(true), 1, function(x) sqrt(sum(x ^ 2))) / sqrt(2)

}

#' @export
Beta_SSE = function(Beta_hat, Beta) {

  sum(((Beta - rowMeans(Beta)) - (Beta_hat - rowMeans(Beta_hat))) ^ 2)

}

#' @export
Beta_FPR = function(Beta_hat, Beta) {

  sum(rowSums(abs(Beta_hat)) > 0 & rowSums(abs(Beta)) == 0) / sum(rowSums(abs(Beta)) == 0)

}

#' @export
Beta_TPR = function(Beta_hat, Beta) {

  sum(rowSums(abs(Beta_hat)) > 0 & rowSums(abs(Beta)) > 0) / sum(rowSums(abs(Beta)) > 0)

}
