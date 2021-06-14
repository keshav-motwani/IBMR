# check that fitted Beta of IBMR_no_Gamma is same as group lasso regression from glmnet

library(IBMR)
library(glmnet)

set.seed(11, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")

number_of_levels = 2
number_per_split = 2
label_levels_per_dataset = list(c(2, 2), c(2, 2), c(2, 2), c(2, 2))
category_mappings = simulate_category_mappings(number_of_levels, number_per_split, label_levels_per_dataset)

p = 50
nonzero = 10

alpha = simulate_alpha(category_mappings$categories, 0.1, 0.5)
Beta = simulate_Beta(category_mappings$categories, p, nonzero, -0.5, 0.5)

X_list = simulate_X_star_list(rep(10000, length(label_levels_per_dataset)), p)
Z_list = simulate_X_star_list(rep(10000, length(label_levels_per_dataset)), p)
Y_list = simulate_Y_list(category_mappings$categories, category_mappings$inverse_category_mappings, X_list, alpha, Beta)

X_list_val = simulate_X_star_list(rep(10000, length(label_levels_per_dataset)), p)
Y_list_val = simulate_Y_list(category_mappings$categories, category_mappings$inverse_category_mappings, X_list_val, alpha, Beta)

system.time({fit = glmnet(do.call(rbind, X_list), unlist(Y_list), family = "multinomial", alpha = 1, standardize = TRUE, intercept = TRUE, type.multinomial = "grouped", nlambda = 20, lambda.min.ratio = 1e-3)})

test = IBMR_no_Gamma(Y_list, category_mappings$categories, category_mappings$category_mappings, X_list, n_lambda = 20, lambda_min_ratio = 1e-3)

plot(coef(fit, fit$lambda[10])[[1]][-1], test$model_fits[[10]]$Beta[, 1])
abline(0, 1)

val = compute_tuning_performance_no_Gamma(test, Y_list_val, category_mappings$categories, category_mappings$category_mappings, X_list_val)

which_min(val)

test = IBMR(Y_list, category_mappings$categories, category_mappings$category_mappings, X_list, X_list, Gamma_update = "gradient")

val = compute_tuning_performance(test, Y_list_val, category_mappings$categories, category_mappings$category_mappings, X_list_val)

which_min(val)
