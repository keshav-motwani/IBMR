# check that objective function decreases

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

X_list = simulate_X_list(rep(10000, length(label_levels_per_dataset)), p)
Y_list = simulate_Y_list(category_mappings$categories, category_mappings$inverse_category_mappings, X_list, alpha, Beta)

Y_list_fine = lapply(Y_list, names)

Y_matrix_list = lapply(1:length(Y_list), function(i) create_Y_matrix(Y_list[[i]], category_mappings$categories, category_mappings$category_mappings[[i]]))
fine_category_mapping = category_mappings$categories
names(fine_category_mapping) = fine_category_mapping
Y_matrix_list_fine = lapply(1:length(Y_list_fine), function(i) create_Y_matrix(Y_list_fine[[i]], category_mappings$categories, as.list(fine_category_mapping)))

system.time({test = fit_alpha_Beta_Gamma(Y_matrix_list, X_list, X_list, 0.01, 0.5, 1000, 1e-6, rep(0, 4), matrix(0, nrow = 20, ncol = 4), lapply(1:length(X_list), function(x) matrix(0, nrow = 20, ncol = 4)))})
system.time({test_fine = fit_alpha_Beta_Gamma(Y_matrix_list_fine, X_list, X_list, 0.01, 0.5, 1000, 1e-6, rep(0, 4), matrix(0, nrow = 20, ncol = 4), lapply(1:length(X_list), function(x) matrix(0, nrow = 20, ncol = 4)))})

all(diff(test$objective[test$objective != 0]) <= 0)
all(diff(test_fine$objective[test_fine$objective != 0]) <= 0)

test = test$Beta
test_fine = test_fine$Beta

par(mfrow = c(2, 2))

plot(test_fine - rowMeans(test_fine), test - rowMeans(test))
abline(0, 1)

plot(test_fine - rowMeans(test_fine), Beta - rowMeans(Beta))
abline(0, 1)

plot(test - rowMeans(test), Beta - rowMeans(Beta))
abline(0, 1)
