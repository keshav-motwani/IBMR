print("Check IBMR and glmnet functions")

TOLERANCE = 1e-12
COEF_THRESHOLD = 1e-4

set.seed(1)

library(IBMR)
library(glmnet)

number_of_levels = 2
splits_per_level = 2
label_levels_per_dataset = list(c(1, 2), c(2, 1), c(1, 2), c(2, 2))
category_mappings = simulate_category_mappings(number_of_levels, splits_per_level, label_levels_per_dataset)

p = 50
nonzero = 10

alpha = simulate_alpha(category_mappings$categories, 0.1, 0.5)
Beta = simulate_Beta(category_mappings$categories, p, nonzero, -0.5, 0.5)

X_list = simulate_X_star_list(rep(100, length(label_levels_per_dataset)), p)
Y_list = simulate_Y_list(category_mappings$categories, category_mappings$inverse_category_mappings, X_list, alpha, Beta)

X_list_val = simulate_X_star_list(rep(100, length(label_levels_per_dataset)), p)
Y_list_val = simulate_Y_list(category_mappings$categories, category_mappings$inverse_category_mappings, X_list_val, alpha, Beta)

system.time({test = IBMR_no_Gamma_subset(Y_list, category_mappings$categories, category_mappings$category_mappings, X_list, Y_list_validation = Y_list_val, category_mappings_validation = category_mappings$category_mappings, X_list_validation = X_list_val, tolerance = TOLERANCE)})
system.time({test2 = glmnet_subset(Y_list, category_mappings$categories, category_mappings$category_mappings, X_list, n_iter = 1e6, Y_list_validation = Y_list_val, category_mappings_validation = category_mappings$category_mappings, X_list_validation = X_list_val, tolerance = TOLERANCE)})

test_that("Estimated Beta from IBMR_no_Gamma_subset matches glmnet_subset", {
  expect(all(abs(test$model_fits[[5]]$Beta[, 1] - test2$model_fits[[5]]$Beta[, 1]) < COEF_THRESHOLD), "coefficients not equal")
})

plot(test$model_fits[[5]]$Beta[, 1], test2$model_fits[[5]]$Beta[, 1])
abline(0, 1)

system.time({test = IBMR_no_Gamma_relabel(Y_list, category_mappings$categories, category_mappings$category_mappings, X_list, Y_list_validation = Y_list_val, category_mappings_validation = category_mappings$category_mappings, X_list_validation = X_list_val, tolerance = TOLERANCE)})
system.time({test2 = glmnet_relabel(Y_list, category_mappings$categories, category_mappings$category_mappings, X_list, Y_list_validation = Y_list_val, category_mappings_validation = category_mappings$category_mappings, X_list_validation = X_list_val, n_iter = 1e6, tolerance = TOLERANCE)})

test_that("Estimated Beta from IBMR_no_Gamma_relabel matches glmnet_relabel", {
  expect(all(abs(test$model_fits[[5]]$Beta[, 1] - test2$model_fits[[5]]$Beta[, 1]) < COEF_THRESHOLD), "coefficients not equal")
})

plot(test$model_fits[[5]]$Beta[, 1], test2$model_fits[[5]]$Beta[, 1])
abline(0, 1)
