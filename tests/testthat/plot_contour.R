library(IBMR)

set.seed(1)

number_of_levels = 2
number_per_split = 2
label_levels_per_dataset = list(c(2, 1), c(2, 2))
category_mappings = simulate_category_mappings(number_of_levels, number_per_split, label_levels_per_dataset)

p = 1
nonzero = 1

alpha = rep(0, 4)
Beta = simulate_Beta(category_mappings$categories, p, nonzero, -0.5, 0.5)

n = c(10, 1000)

X_list = simulate_X_star_list(n, p)
Y_list = simulate_Y_list(category_mappings$categories, category_mappings$inverse_category_mappings, X_list, alpha, Beta)
Y_list_fine = lapply(Y_list, names)

Y_matrix_list = lapply(1:length(Y_list), function(i) create_Y_matrix(Y_list[[i]], category_mappings$categories, category_mappings$category_mappings[[i]]))

fine_category_mapping = category_mappings$categories
names(fine_category_mapping) = fine_category_mapping
Y_matrix_list_fine = lapply(1:length(Y_list_fine), function(i) create_Y_matrix(Y_list_fine[[i]], category_mappings$categories, as.list(fine_category_mapping)))

range = 5

Beta_grid = expand.grid(Beta_1 = seq(-range, range, length.out = 100), Beta_2 = seq(-range, range, length.out = 100))
Beta_grid$Beta_3 = Beta[3]
Beta_grid$Beta_4 = Beta[4]

nll = apply(Beta_grid, 1, function(x) IBMR:::compute_negative_log_likelihood(Y_matrix_list, X_list, X_list, alpha, matrix(x, nrow = 1), lapply(1:2, function(x) matrix(0, nrow = 1, ncol = 4)), sum(n)))
# nll_fine = apply(Beta_grid, 1, function(x) IBMR:::compute_negative_log_likelihood(Y_matrix_list_fine, X_list, X_list, alpha, matrix(x, nrow = 1), lapply(1:2, function(x) matrix(0, nrow = 1, ncol = 4)), sum(n)))

filled.contour(x = seq(-range, range, length.out = 100), y = seq(-range, range, length.out = 100), z = matrix(nll, nrow = 100, byrow = TRUE))
# filled.contour(x = seq(-range, range, length.out = 100), y = seq(-range, range, length.out = 100), z = matrix(nll_fine, nrow = 100, byrow = TRUE))
