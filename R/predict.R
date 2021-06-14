#' @export
predict_categories = function(X, alpha, Beta, categories) {

  P = compute_probabilities_no_Gamma(X, alpha, Beta)

  categories = categories[apply(P, 1, which.max)]

  return(categories)

}
