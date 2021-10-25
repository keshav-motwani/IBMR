#' @export
simulate_category_mappings = function(number_of_levels, splits_per_level, label_levels_per_dataset) {

  # if(!all(sapply(1:splits_per_level, function(i) number_of_levels %in% sapply(label_levels_per_dataset, `[`, i)))) {
  #   stop("The finest level categories (highest level) must be present in at least one dataset as specified in label_levels_per_dataset.")
  # }

  if (length(splits_per_level) == 1) splits_per_level = rep(splits_per_level, number_of_levels)

  if (number_of_levels != length(splits_per_level)) stop("length of splits_per_level must be the same as number_of_levels")

  if (!all(sapply(label_levels_per_dataset, length) == splits_per_level[1])) {
    stop("For each dataset, one level for each of the first branches of the tree must be specified (for a total of splits_per_level categories).")
  }

  label_tree = create_label_tree(number_of_levels, splits_per_level)

  category_mappings = lapply(label_levels_per_dataset, function(x) get_category_mapping(label_tree, x))

  return(list(category_mappings = lapply(category_mappings, `[[`, "category_mapping"),
              inverse_category_mappings = lapply(category_mappings, `[[`, "inverse_category_mapping"),
              categories = label_tree[nrow(label_tree), ],
              number_of_levels = number_of_levels,
              splits_per_level = splits_per_level))

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

create_label_tree = function(number_of_levels, splits_per_level) {

  number_of_categories = prod(splits_per_level)

  tree = matrix(nrow = number_of_levels, ncol = number_of_categories)

  tree[1, ] = rep(1:splits_per_level[1], each = number_of_categories / splits_per_level[1])

  for (i in 2:number_of_levels) {
    tree[i, ] = rep(rep(1:splits_per_level[i], prod(splits_per_level[1:(i - 1)])), each = number_of_categories / prod(splits_per_level[1:i]))
  }

  label_tree = t(sapply(2:number_of_levels, function(i) apply(tree[1:i, ], 2, function(x) paste(x, collapse = ""))))
  label_tree = rbind(tree[1, ], label_tree)

  return(label_tree)

}

#' @export
simulate_X_star = function(n, p, mean = 0, rho = 0.5) {

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
simulate_X_star_list = function(n_k, p, mean = 0, rho = 0.5) {

  X = simulate_X_star(sum(n_k), p, mean, rho)

  lapply(split(X, rep(1:length(n_k), n_k)), matrix, ncol = p)

}

#' @export
simulate_U = function(X_star, rank, batch_effect) {

  n = nrow(X_star)
  p = ncol(X_star)

  if (rank == "int") {

    U = matrix(rep(1, nrow(X_star)), ncol = 1) %*% matrix(rnorm(p), nrow = 1)

  } else {

    rank = as.numeric(rank)
    U = matrix(rnorm(n * rank), ncol = rank) %*% matrix(rnorm(rank * p), nrow = rank)

  }

  c = batch_effect * norm(X_star, "F") / norm(U, "F")

  U = c * U

  print(norm(U, "F") / norm(X_star, "F"))

  return(U)

}

#' @export
simulate_U_list = function(X_star_list, rank, batch_effect) {

  U_list = lapply(X_star_list, simulate_U, rank = rank, batch_effect = batch_effect)

  return(U_list)

}

#' @export
compute_X_list = function(X_star_list, U_list) {

  mapply(`+`, X_star_list, U_list, SIMPLIFY = FALSE)

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
simulate_structured_Beta = function(splits_per_level, p, nonzero, features_per_level, sigma = 2) {

  stopifnot(sum(features_per_level) == nonzero)

  indices = sample(1:p, nonzero)

  K = prod(splits_per_level)
  Beta = matrix(0, nrow = p, ncol = K)

  for (l in 1:length(splits_per_level)) {

    level_indices = sample(indices, features_per_level[l])
    print(level_indices)

    for (c in 1:prod(splits_per_level[1:l])) {

      Beta[level_indices, (((c - 1) * K / prod(splits_per_level[1:l])) + 1):(c * K / prod(splits_per_level[1:l]))] = rnorm(length(level_indices), sd = sigma) %*% t(rep(1, K / prod(splits_per_level[1:l])))

    }

    indices = setdiff(indices, level_indices)

  }

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
