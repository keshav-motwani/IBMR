#' @export
simulate_category_mappings = function(number_of_levels, number_per_split, label_levels_per_dataset) {

  # if(!all(sapply(1:number_per_split, function(i) number_of_levels %in% sapply(label_levels_per_dataset, `[`, i)))) {
  #   stop("The finest level categories (highest level) must be present in at least one dataset as specified in label_levels_per_dataset.")
  # }

  if (!all(sapply(label_levels_per_dataset, length) == number_per_split)) {
    stop("For each dataset, one level for each of the first branches of the tree must be specified (for a total of number_per_split categories).")
  }

  number_of_categories = number_per_split ^ number_of_levels

  tree = matrix(nrow = number_of_levels, ncol = number_of_categories)

  for (i in 1:number_of_levels) {
    tree[i, ] = rep(rep(1:number_per_split, number_per_split ^ (i - 1)), each = number_per_split ^ (number_of_levels - i))
  }

  label_tree = t(sapply(2:number_of_levels, function(i) apply(tree[1:i, ], 2, function(x) paste(x, collapse = ""))))
  label_tree = rbind(tree[1, ], label_tree)

  category_mappings = lapply(label_levels_per_dataset, function(x) get_category_mapping(label_tree, x))

  return(list(category_mappings = lapply(category_mappings, `[[`, "category_mapping"),
              inverse_category_mappings = lapply(category_mappings, `[[`, "inverse_category_mapping"),
              categories = label_tree[nrow(label_tree), ]))

}

get_category_mapping = function(label_tree, label_levels) {

  inverse_category_mapping = c()

  for (label in 1:length(label_levels)) {

    label_mapping = label_tree[label_levels[label], label_tree[1, ] == as.character(label)]
    names(label_mapping) = label_tree[nrow(label_tree), label_tree[1, ] == as.character(label)]

    inverse_category_mapping = c(inverse_category_mapping, label_mapping)

  }

  category_mapping = list()

  for (label in unique(inverse_category_mapping)) {

    category_mapping[[label]] = names(inverse_category_mapping)[which(inverse_category_mapping == label)]

  }

  return(list(category_mapping = category_mapping, inverse_category_mapping = inverse_category_mapping))

}

#' @export
simulate_X = function(n, p, mean = 0, rho = 0.5) {

  Sigma = matrix(nrow = p, ncol = p)

  for (i in 1:p) {
    for (j in 1:p) {
      Sigma[i, j] = rho ^ abs(i - j)
    }
  }

  eo = eigen(Sigma)

  X = matrix(rnorm(n * nrow(Sigma)), nrow = n) %*% eo$vec %*% diag(eo$val^.5) %*% t(eo$vec)

  X = X + mean

  return(X)

}

#' @export
simulate_X_list = function(n_k, p, mean = 0, rho = 0.5) {

  lapply(n_k, simulate_X, p, mean, rho)

}

#' @export
simulate_alpha = function(categories, lower = -2, upper = 2) {

  alpha = runif(length(categories), lower, upper)

  names(alpha) = categories

  return(alpha - mean(alpha))

}

#' @export
simulate_Beta = function(categories, p, nonzero, lower = -2, upper = 2) {

  Beta = matrix(0, nrow = p, ncol = length(categories))

  Beta[1:nonzero, ] = runif(nonzero * length(categories), lower, upper)

  colnames(Beta) = categories

  return(Beta)

}

#' @export
simulate_Y = function(categories, inverse_category_mapping, X, alpha, Beta) {

  P = compute_probabilities_no_Gamma(X, alpha, Beta)

  Y_fine = apply(P, 1, function(x) sample(categories, 1, prob = x))
  Y = inverse_category_mapping[Y_fine]

  return(Y)

}

#' @export
simulate_Y_list = function(categories, inverse_category_mappings, X_list, alpha, Beta) {

  lapply(1:length(X_list), function(i) simulate_Y(categories, inverse_category_mappings[[i]], X_list[[i]], alpha, Beta))

}
