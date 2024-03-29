print("Check that fit_alpha_Beta_Gamma satisfies KKT conditions and objective function always decreases with coarse data")

TOLERANCE = 1e-12
KKT_THRESHOLD = 1e-6
COEF_THRESHOLD = 1e-4

set.seed(1)

library(IBMR)
library(glmnet)
library(CVXR)

number_of_levels = 2
splits_per_level = 2
label_levels_per_dataset = list(c(1, 2), c(2, 1))
category_mappings = simulate_category_mappings(number_of_levels, splits_per_level, label_levels_per_dataset)

p = 20
nonzero = 10

alpha = simulate_alpha(category_mappings$categories, 0.1, 0.5)
Beta = simulate_Beta(category_mappings$categories, p, nonzero, -0.5, 0.5)

X_list = simulate_X_star_list(rep(100, length(label_levels_per_dataset)), p)
Z_list = simulate_X_star_list(rep(100, length(label_levels_per_dataset)), p)
Y_list = simulate_Y_list(category_mappings$categories, category_mappings$inverse_category_mappings, X_list, alpha, Beta)

Y_matrix_list = lapply(1:length(Y_list), function(i) create_Y_matrix(Y_list[[i]], category_mappings$categories, category_mappings$category_mappings[[i]]))

lambda = 0.01
rho = 0.5

system.time({test = fit_alpha_Beta_Gamma(Y_matrix_list, X_list, Z_list, lambda, rho, 1000, TOLERANCE, rep(0, 4), matrix(0, nrow = 20, ncol = 4), lapply(1:length(X_list), function(x) matrix(0, nrow = 20, ncol = 4)))})
system.time({test_Newton = fit_alpha_Beta_Gamma_Newton(Y_matrix_list, X_list, Z_list, lambda, rho, 1000, TOLERANCE, rep(0, 4), matrix(0, nrow = 20, ncol = 4), lapply(1:length(X_list), function(x) matrix(0, nrow = 20, ncol = 4)))})

test_that("Estimate from fit_alpha_Beta_Gamma for coarse resolution data satisfies KKT conditions (necessary, but not sufficient for optimality as nonconvex)", {
  expect(all(check_KKT_IBMR(Y_matrix_list, X_list, Z_list, lambda, rho, test$alpha, test$Beta, test$Gamma_list) - c(0, 0, 1, 0) < KKT_THRESHOLD), "doesn't satisfy KKT")
})

test_that("Estimate from fit_alpha_Beta_Gamma_Newton for coarse resolution data satisfies KKT conditions (necessary, but not sufficient for optimality as nonconvex)", {
  expect(all(check_KKT_IBMR(Y_matrix_list, X_list, Z_list, lambda, rho, test_Newton$alpha, test_Newton$Beta, test_Newton$Gamma_list) - c(0, 0, 1, 0) < KKT_THRESHOLD), "doesn't satisfy KKT")
})

test_that("Objective function for fit_alpha_Beta_Gamma is always decreasing", {
  expect(all(diff(test$objective[test$objective != 0]) <= 0), "objective function increased somewhere")
})

test_that("Objective function for fit_alpha_Beta_Gamma_Newton is always decreasing", {
  expect(all(diff(test_Newton$objective[test_Newton$objective != 0]) <= 0), "objective function increased somewhere")
})
