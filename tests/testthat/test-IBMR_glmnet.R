print("Check IBMR and glmnet functions")

TOLERANCE = 1e-12
COEF_THRESHOLD = 1e-4

set.seed(1)

library(IBMR)
library(glmnet)

number_of_levels = 2
number_per_split = 2
label_levels_per_dataset = list(c(2, 2), c(2, 2), c(2, 2), c(2, 2))
category_mappings = simulate_category_mappings(number_of_levels, number_per_split, label_levels_per_dataset)

p = 50
nonzero = 10

alpha = simulate_alpha(category_mappings$categories, 0.1, 0.5)
Beta = simulate_Beta(category_mappings$categories, p, nonzero, -0.5, 0.5)

X_list = simulate_X_star_list(rep(100, length(label_levels_per_dataset)), p)
Z_list = simulate_X_star_list(rep(100, length(label_levels_per_dataset)), p)
Y_list = simulate_Y_list(category_mappings$categories, category_mappings$inverse_category_mappings, X_list, alpha, Beta)

X_list_val = simulate_X_star_list(rep(100, length(label_levels_per_dataset)), p)
Y_list_val = simulate_Y_list(category_mappings$categories, category_mappings$inverse_category_mappings, X_list_val, alpha, Beta)

fit = glmnet(do.call(rbind, X_list), unlist(Y_list), family = "multinomial", alpha = 1, standardize = TRUE, intercept = TRUE, type.multinomial = "grouped", nlambda = 25, lambda.min.ratio = 1e-4, thresh = TOLERANCE)

system.time({test = IBMR_no_Gamma(Y_list, category_mappings$categories, category_mappings$category_mappings, X_list, tolerance = TOLERANCE)})

test_that("Estimated Beta from IBMR_no_Gamma matches glmnet for fine resolution data", {
  expect(all(abs(coef(fit, fit$lambda[10])[[1]][-1] - test$model_fits[[10]]$Beta[, 1]) < COEF_THRESHOLD), "coefficients not equal")
})

plot(coef(fit, fit$lambda[10])[[1]][-1], test$model_fits[[10]]$Beta[, 1])
abline(0, 1)

# should be equal to IBMR_no_Gamma when all are at finest resolution
test2 = glmnet_subset(Y_list, category_mappings$categories, category_mappings$category_mappings, X_list, n_lambda = 25, lambda_min_ratio = 1e-4, tolerance = TOLERANCE)

test_that("Estimated Beta from glmnet_subset matches glmnet for fine resolution data", {
  expect(all(abs(coef(fit, fit$lambda[10])[[1]][-1] - test2$model_fits[[10]]$Beta[, 1]) < COEF_THRESHOLD), "coefficients not equal")
})

plot(coef(fit, fit$lambda[10])[[1]][-1], test2$model_fits[[10]]$Beta[, 1])
abline(0, 1)

# should be equal to IBMR_no_Gamma when all are at finest resolution
test = glmnet_relabel(Y_list, category_mappings$categories, category_mappings$category_mappings, X_list, n_rho = 20, rho_min_ratio = 1e-3, n_lambda = 25, lambda_min_ratio = 1e-4, tolerance = TOLERANCE)

test_that("Estimated Beta from glmnet_relabel matches glmnet for fine resolution data", {
  expect(all(abs(coef(fit, fit$lambda[10])[[1]][-1] - test$model_fits[[4]][[10]]$Beta[, 1]) < COEF_THRESHOLD), "coefficients not equal")
})

plot(coef(fit, fit$lambda[10])[[1]][-1], test$model_fits[[4]][[10]]$Beta[, 1])
abline(0, 1)

fit = cv.glmnet(do.call(rbind, X_list), unlist(Y_list), family = "multinomial", alpha = 1, standardize = TRUE, intercept = TRUE, type.multinomial = "grouped", nlambda = 25, lambda.min.ratio = 1e-4)

test = glmnet_split(list(unlist(Y_list)), category_mappings$categories, category_mappings$category_mappings, list(do.call(rbind, X_list)), n_lambda = 25, lambda_min_ratio = 1e-4)

probs_cv_glmnet = predict(fit, X_list[[1]], type = "response", s = fit$lambda.min)[, , 1]

probs_glmnet_split = predict_probabilities_glmnet_split(test, list(X_list[[1]]))[[1]]

test_that("Estimated probabilities from glmnet_split matches cv.glmnet for fine resolution data with 1 dataset", {
  expect(all(abs(probs_cv_glmnet - probs_glmnet_split) < COEF_THRESHOLD), "probabilities not equal")
})

fit = cv.glmnet(X_list[[1]], Y_list[[1]], family = "multinomial", alpha = 1, standardize = TRUE, intercept = TRUE, type.multinomial = "grouped", nlambda = 25, lambda.min.ratio = 1e-4)
fit2 = cv.glmnet(X_list[[2]], Y_list[[2]], family = "multinomial", alpha = 1, standardize = TRUE, intercept = TRUE, type.multinomial = "grouped", nlambda = 25, lambda.min.ratio = 1e-4)

test = glmnet_split(Y_list[1:2], category_mappings$categories, category_mappings$category_mappings[1:2], X_list[1:2], n_lambda = 25, lambda_min_ratio = 1e-4)

probs_cv_glmnet = predict(fit, X_list[[1]], type = "response", s = fit$lambda.min)[, , 1]
probs_cv_glmnet2 = predict(fit2, X_list[[1]], type = "response", s = fit2$lambda.min)[, , 1]
probs_cv_glmnet = (probs_cv_glmnet + probs_cv_glmnet2) / 2

probs_glmnet_split = predict_probabilities_glmnet_split(test, list(X_list[[1]]))[[1]]

test_that("Estimated probabilities from glmnet_split matches cv.glmnet for fine resolution data with 2 datasets", {
  expect(all(abs(probs_cv_glmnet - probs_glmnet_split) < COEF_THRESHOLD), "probabilities not equal")
})

plot(probs_cv_glmnet, probs_glmnet_split)
abline(0, 1)

system.time({test_no_Gamma = IBMR_no_Gamma(Y_list, category_mappings$categories, category_mappings$category_mappings, X_list, tolerance = TOLERANCE, n_lambda = 10)})
system.time({test = IBMR(Y_list, category_mappings$categories, category_mappings$category_mappings, X_list, X_list, Gamma_update = "gradient", n_lambda = 10, n_rho = 5, rho_min_ratio = 1e-4)})

plot(test_no_Gamma$model_fits[[10]]$Beta, test$model_fits[[1]][[10]]$Beta)
abline(0, 1)

test_that("Estimated Beta from IBMR approx matches IBMR_no_Gamma for largest rho", {
  expect(max(abs(test_no_Gamma$model_fits[[10]]$Beta - test$model_fits[[1]][[10]]$Beta)) < 0.005, "coefficients not equal")
})
