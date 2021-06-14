# check that fitted Beta of fit_alpha_Beta is same as group lasso regression from glmnet when Gamma = 0

library(IBMR)
library(glmnet)

set.seed(11, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")

number_of_levels = 2
number_per_split = 2
label_levels_per_dataset = list(c(1, 1), c(1, 2), c(2, 1), c(1, 1))
category_mappings = simulate_category_mappings(number_of_levels, number_per_split, label_levels_per_dataset)

p = 20
nonzero = 20

alpha = simulate_alpha(category_mappings$categories, 0.1, 0.5)
Beta = simulate_Beta(category_mappings$categories, p, nonzero, -0.5, 0.5)

X_list = simulate_X_star_list(rep(10000, length(label_levels_per_dataset)), p)
Y_list = simulate_Y_list(category_mappings$categories, category_mappings$inverse_category_mappings, X_list, alpha, Beta)

Y_list_fine = lapply(Y_list, names)

Y_matrix_list = lapply(1:length(Y_list), function(i) create_Y_matrix(Y_list[[i]], category_mappings$categories, category_mappings$category_mappings[[i]]))
fine_category_mapping = category_mappings$categories
names(fine_category_mapping) = fine_category_mapping
Y_matrix_list_fine = lapply(1:length(Y_list_fine), function(i) create_Y_matrix(Y_list_fine[[i]], category_mappings$categories, as.list(fine_category_mapping)))

system.time({fit3 = glmnet(do.call(rbind, X_list), unlist(Y_list_fine), family = "multinomial", alpha = 1, standardize = FALSE, intercept = TRUE, type.multinomial = "grouped")})
test2 = as.matrix(do.call(cbind, coef(fit3, fit3$lambda[20]))[-1, ])

system.time({test = fit_alpha_Beta(Y_matrix_list, X_list, X_list, fit3$lambda[20], 1000, 1e-6, rep(0, 4), matrix(0, nrow = 20, ncol = 4))$Beta})
system.time({test_fine = fit_alpha_Beta(Y_matrix_list_fine, X_list, X_list, fit3$lambda[20], 1000, 1e-6, rep(0, 4), matrix(0, nrow = 20, ncol = 4))})

fit3$a0[, 20]
test_fine$alpha

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

test_fine = fit_alpha_Beta(Y_matrix_list_fine, X_list, X_list, fit3$lambda[1], 10000, 1e-12, rep(0, 4), matrix(0, nrow = 20, ncol = 4))
fit3$a0[, 1]
test_fine$alpha
test_fine$Beta

test_fine = fit_alpha_Beta(Y_matrix_list_fine, X_list, X_list, fit3$lambda[1] * 1.000001, 10000, 1e-12, rep(0, 4), matrix(0, nrow = 20, ncol = 4))
fit3$a0[, 1]
test_fine$alpha
test_fine$Beta

rm(list = ls())
