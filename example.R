library(IBMR)

set.seed(1)

p = 500
nonzero = 100
n = rep(500, 2)

category_mappings = simulate_category_mappings(2, 2, list(c(1, 2), c(2, 1)))

alpha = simulate_alpha(category_mappings$categories)
Beta = simulate_Beta(category_mappings$categories, p, nonzero)

X_star_list = simulate_X_star_list(n, p)
U_list = simulate_U_list(X_star_list, "int", 0.1)
X_list = compute_X_list(X_star_list, U_list)
Y_list = simulate_Y_list(category_mappings$categories, category_mappings$inverse_category_mappings, X_star_list, alpha, Beta)
Z_list = list(matrix(1, nrow = n[1]), matrix(1, nrow = n[2]))

X_star_list_val = simulate_X_star_list(n, p)
U_list_val = simulate_U_list(X_star_list_val, "int", 0.1)
X_list_val = compute_X_list(X_star_list_val, U_list_val)
Y_list_val = simulate_Y_list(category_mappings$categories, category_mappings$inverse_category_mappings, X_star_list_val, alpha, Beta)

X_list_test = simulate_X_star_list(10000, p)
Y_list_test = simulate_Y_list(category_mappings$categories, list(setNames(nm = category_mappings$categories)), X_list_test, alpha, Beta)

system.time({IBMR_fit = IBMR(Y_list = Y_list,
                             categories = category_mappings$categories,
                             category_mappings = category_mappings$category_mappings,
                             X_list = X_list,
                             Z_list = Z_list,
                             Y_list_validation = Y_list_val,
                             category_mappings_validation = category_mappings$category_mappings,
                             X_list_validation = X_list_val,
                             verbose = FALSE)})

system.time({IBMR_NG_fit = IBMR_no_Gamma(Y_list = Y_list,
                                         categories = category_mappings$categories,
                                         category_mappings = category_mappings$category_mappings,
                                         X_list = X_list,
                                         Y_list_validation = Y_list_val,
                                         category_mappings_validation = category_mappings$category_mappings,
                                         X_list_validation = X_list_val,
                                         verbose = FALSE)})

system.time({subset_fit = subset(Y_list = Y_list,
                                 categories = category_mappings$categories,
                                 category_mappings = category_mappings$category_mappings,
                                 X_list = X_list,
                                 Y_list_validation = Y_list_val,
                                 category_mappings_validation = category_mappings$category_mappings,
                                 X_list_validation = X_list_val,
                                 verbose = FALSE)})

system.time({relabel_fit = relabel(Y_list = Y_list,
                                   categories = category_mappings$categories,
                                   category_mappings = category_mappings$category_mappings,
                                   X_list = X_list,
                                   Y_list_validation = Y_list_val,
                                   category_mappings_validation = category_mappings$category_mappings,
                                   X_list_validation = X_list_val,
                                   verbose = FALSE)})

mean(Y_list_test[[1]] != predict_categories(predict_probabilities(IBMR_fit$best_model, X_list_test))[[1]])
mean(Y_list_test[[1]] != predict_categories(predict_probabilities(IBMR_NG_fit$best_model, X_list_test))[[1]])
mean(Y_list_test[[1]] != predict_categories(predict_probabilities(subset_fit$best_model, X_list_test))[[1]])
mean(Y_list_test[[1]] != predict_categories(predict_probabilities(relabel_fit$best_model, X_list_test))[[1]])

