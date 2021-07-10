print("Check that probabilities from fit_alpha match sample proportions when all data is at finest resolution")

TOLERANCE = 1e-12
PROB_THRESHOLD = 1e-4

set.seed(1)

library(IBMR)
library(glmnet)
library(CVXR)

number_of_levels = 2
number_per_split = 2
label_levels_per_dataset = list(c(2, 2), c(2, 2))
category_mappings = simulate_category_mappings(number_of_levels, number_per_split, label_levels_per_dataset)

p = 20
nonzero = 10

alpha = simulate_alpha(category_mappings$categories, 0.1, 0.5)
Beta = simulate_Beta(category_mappings$categories, p, nonzero, -0.5, 0.5)

X_list = simulate_X_star_list(rep(100, length(label_levels_per_dataset)), p)
Y_list = simulate_Y_list(category_mappings$categories, category_mappings$inverse_category_mappings, X_list, alpha, Beta)

Y_matrix_list = lapply(1:length(Y_list), function(i) create_Y_matrix(Y_list[[i]], category_mappings$categories, category_mappings$category_mappings[[i]]))

system.time({test = fit_alpha(Y_matrix_list, X_list, X_list, 1000, TOLERANCE, rep(0, 4))})

test_that("Probabilities from fit_alpha for fine resolution data match sample proportions", {
  expect(all(abs(exp(test$alpha)/sum(exp(test$alpha)) - colMeans(do.call(rbind, Y_matrix_list))) < PROB_THRESHOLD), "probabilities don't match")
})

test_that("Objective function for fit_alpha is always decreasing", {
  expect(all(diff(test$objective[test$objective != 0]) <= 1e-12), "objective function increased somewhere")
})
