#' @export
create_Y_matrix = function(Y, categories, category_mapping) {

  Y = t(sapply(Y, function(y) categories %in% category_mapping[[y]]) * 1)

  colnames(Y) = categories

  return(Y)

}
