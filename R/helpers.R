#' @export
subset_helper = function(Y_list,
                         categories,
                         category_mappings,
                         X_list,
                         combine = TRUE) {

  Y_matrix_list = lapply(1:length(Y_list), function(i) create_Y_matrix(Y_list[[i]], categories, category_mappings[[i]]))

  indices_list = lapply(Y_matrix_list, function(Y) which(rowSums(Y) == 1))

  print(paste0("Keeping ", sum(sapply(indices_list, length)), " observations out of a total of ", sum(sapply(Y_list, length))))

  Y_list = mapply(Y = Y_list, indices = indices_list, map = category_mappings, FUN = function(Y, indices, map) unlist(map[Y[indices]]), SIMPLIFY = FALSE)
  X_list = mapply(X = X_list, indices = indices_list, FUN = function(X, indices) X[indices, ], SIMPLIFY = FALSE)

  if (combine) {
    Y_list = list(unlist(Y_list))
    X_list = list(do.call(rbind, X_list))
  }

  return(list(Y_list = Y_list, X_list = X_list))

}
