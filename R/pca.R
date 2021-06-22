#' @export
compute_pca_for_Z_list = function(X_list, n_pcs, intercept = TRUE) {

  if (n_pcs > 0) {

    set.seed(1)

    k = rep(1:length(X_list), sapply(X_list, nrow))

    X = do.call(rbind, X_list)
    X = scale(X, scale = FALSE)
    Z = irlba::irlba(X, n_pcs)
    Z = X %*% Z$v

    if (intercept) Z = cbind(1, Z)

    Z_list = lapply(1:length(X_list), function(i) Z[k == i, ])

  } else {

    Z_list = lapply(X_list, function(x) matrix(1, nrow = nrow(x), ncol = 1))

  }

  return(Z_list)

}
