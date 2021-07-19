#' @export
create_Y_matrix = function(Y, categories, category_mapping) {

  Y = t(sapply(Y, function(y) categories %in% category_mapping[[y]]) * 1)

  colnames(Y) = categories

  return(Y)

}

#' @export
create_fine_category_mappings = function(categories, K) {

  names(categories) = categories

  list(category_mappings = replicate(K, as.list(categories), simplify = FALSE),
       inverse_category_mappings = replicate(K, categories, simplify = FALSE),
       categories = categories)

}
