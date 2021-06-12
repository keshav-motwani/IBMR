# check fit_alpha against sample proportions for fine case

library(IBMR)
library(glmnet)
library(CVXR)

set.seed(11, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")

number_of_levels = 2
number_per_split = 2
label_levels_per_dataset = list(c(2, 2), c(2, 2))
category_mappings = simulate_category_mappings(number_of_levels, number_per_split, label_levels_per_dataset)

p = 20
nonzero = 10

alpha = simulate_alpha(category_mappings$categories, 0.1, 0.5)
Beta = simulate_Beta(category_mappings$categories, p, nonzero, -0.5, 0.5)

X_list = simulate_X_list(rep(100, length(label_levels_per_dataset)), p)
Y_list = simulate_Y_list(category_mappings$categories, category_mappings$inverse_category_mappings, X_list, alpha, Beta)

Y_matrix_list = lapply(1:length(Y_list), function(i) create_Y_matrix(Y_list[[i]], category_mappings$categories, category_mappings$category_mappings[[i]]))

system.time({test = fit_alpha(Y_matrix_list, X_list, X_list, 1000, 1e-12, rep(0, 4))})

all(diff(test$objective[test$objective != 0]) <= 0)

all(abs(exp(test$alpha)/sum(exp(test$alpha)) - colMeans(do.call(rbind, Y_matrix_list))) < 1e-8)
