#' @export
subset_helper = function(Y_list,
                         categories,
                         category_mappings,
                         X_list) {

  Y_matrix_list = lapply(1:length(Y_list), function(i) create_Y_matrix(Y_list[[i]], categories, category_mappings[[i]]))

  indices_list = lapply(Y_matrix_list, function(Y) which(rowSums(Y) == 1))

  print(paste0("Keeping ", sum(sapply(indices_list, length)), " observations out of a total of ", sum(sapply(Y_list, length))))

  Y_list = list(unlist(mapply(Y = Y_list, indices = indices_list, map = category_mappings, FUN = function(Y, indices, map) unlist(map[Y[indices]]), SIMPLIFY = FALSE)))
  X_list = list(do.call(rbind, mapply(X = X_list, indices = indices_list, FUN = function(X, indices) X[indices, ], SIMPLIFY = FALSE)))

  return(list(Y_list = Y_list, X_list = X_list))

}
