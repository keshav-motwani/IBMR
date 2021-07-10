print("Check that fit_alpha_Beta matches glmnet group lasso regression estimate when Gamma = 0")

TOLERANCE = 1e-12
KKT_THRESHOLD = 1e-5
COEF_THRESHOLD = 1e-4

set.seed(1)

library(IBMR)
library(glmnet)

number_of_levels = 2
number_per_split = 2
label_levels_per_dataset = list(c(1, 1), c(1, 2), c(2, 1), c(1, 1))
category_mappings = simulate_category_mappings(number_of_levels, number_per_split, label_levels_per_dataset)

p = 20
nonzero = 10

alpha = simulate_alpha(category_mappings$categories)
Beta = simulate_Beta(category_mappings$categories, p, nonzero)

X_list = simulate_X_star_list(rep(200, length(label_levels_per_dataset)), p)
Y_list = simulate_Y_list(category_mappings$categories, category_mappings$inverse_category_mappings, X_list, alpha, Beta)

Y_list_fine = lapply(Y_list, names)

Y_matrix_list = lapply(1:length(Y_list), function(i) create_Y_matrix(Y_list[[i]], category_mappings$categories, category_mappings$category_mappings[[i]]))

fine_category_mapping = category_mappings$categories
names(fine_category_mapping) = fine_category_mapping
Y_matrix_list_fine = lapply(1:length(Y_list_fine), function(i) create_Y_matrix(Y_list_fine[[i]], category_mappings$categories, as.list(fine_category_mapping)))

system.time({fit3 = glmnet(do.call(rbind, X_list), unlist(Y_list_fine), family = "multinomial", alpha = 1, standardize = FALSE, intercept = TRUE, type.multinomial = "grouped", maxit = 1e7, thresh = TOLERANCE)})
test2 = as.matrix(do.call(cbind, coef(fit3, fit3$lambda[p]))[-1, ])

system.time({test = fit_alpha_Beta(Y_matrix_list, X_list, X_list, fit3$lambda[p], 1000, TOLERANCE, rep(0, 4), matrix(0, nrow = p, ncol = 4))})
system.time({test_fine = fit_alpha_Beta(Y_matrix_list_fine, X_list, X_list, fit3$lambda[p], 1000, TOLERANCE, rep(0, 4), matrix(0, nrow = p, ncol = 4))})

test_that("Estimated Beta from fit_alpha_Beta matches glmnet for fine resolution data", {
  expect(all(abs((test_fine$Beta[, 1] - test_fine$Beta[, 4]) - (test2[, 1] - test2[, 4])) < COEF_THRESHOLD), "coefficients not equal")
})

test_that("Estimate from fit_alpha_Beta for fine resolution data satisfies KKT conditions (sufficient for optimality as convex)", {
  expect(all(check_KKT_IBMR_no_Gamma(Y_matrix_list_fine, X_list, fit3$lambda[p], test_fine$alpha, test_fine$Beta) - c(0, 0, 1) < KKT_THRESHOLD), "doesn't satisfy KKT")
})

test_that("Estimate from fit_alpha_Beta for coarse resolution data satisfies KKT conditions (necessary, but not sufficient for optimality as nonconvex)", {
  expect(all(check_KKT_IBMR_no_Gamma(Y_matrix_list, X_list, fit3$lambda[p], test$alpha, test$Beta) - c(0, 0, 1) < KKT_THRESHOLD), "doesn't satisfy KKT")
})

test = test$Beta
test_fine = test_fine$Beta

par(mfrow = c(2, 2))

plot(test_fine - rowMeans(test_fine), test2 - rowMeans(test2))
abline(0, 1)

plot(test_fine - rowMeans(test_fine), test - rowMeans(test))
abline(0, 1)

plot(test_fine - rowMeans(test_fine), Beta - rowMeans(Beta))
abline(0, 1)

plot(test - rowMeans(test), Beta - rowMeans(Beta))
abline(0, 1)

test_fine = fit_alpha_Beta(Y_matrix_list_fine, X_list, X_list, fit3$lambda[1] * 1.000001, 10000, TOLERANCE, rep(0, 4), matrix(0, nrow = p, ncol = 4))
test_fine$Beta

test_that("Estimate of Beta from fit_alpha_Beta for large lambda is 0", {
  expect(all(test_fine$Beta == 0), "nonzero estimate of Beta")
})
