print("Check that fit_Gamma matches glmnet ridge regression estimate when only 1 dataset")

set.seed(1)

library(IBMR)
library(glmnet)

create_Y_matrix_fine = function(Y, categories) {

  Y = t(sapply(Y, function(y) categories == y) * 1)

  colnames(Y) = categories

  return(Y)

}

n = 10000

g4 = sample(1:4, n, replace = TRUE)
x = matrix(rnorm(n * 100), n, 100)

system.time({glmnet_fit = glmnet(x, g4, family = "multinomial", alpha = 0, standardize = FALSE, intercept = FALSE, lambda = 0.5)})
test2 = as.matrix(do.call(cbind, coef(glmnet_fit, 0))[-1, ])

Y_mat = create_Y_matrix_fine(g4, 1:4)

system.time({test = fit_Gamma(list(Y_mat), list(x), list(x), 0.5, 1000, 1e-6, list(matrix(0, ncol = 4, nrow = 100)))$Gamma_list[[1]]})
plot(test[, 1] - test[, 4], test2[, 1] - test2[, 4])
abline(0, 1)

test_that("Estimated Gamma_1 from fit_Gamma matches glmnet", {
  expect(all(abs((test[, 1] - test[, 4]) - (test2[, 1] - test2[, 4])) < 1e-4), "coefficients not equal")
})

system.time({test = fit_Gamma_Newton(list(Y_mat), list(x), list(x), 0.5, 1000, 1e-6, list(matrix(0, ncol = 4, nrow = 100)))$Gamma_list[[1]]})
plot(test[, 1] - test[, 4], test2[, 1] - test2[, 4])
abline(0, 1)

test_that("Estimated Gamma_1 from fit_Gamma_Newton matches glmnet", {
  expect(all(abs((test[, 1] - test[, 4]) - (test2[, 1] - test2[, 4])) < 1e-4), "coefficients not equal")
})

# microbenchmark::microbenchmark(glmnet = {glmnet_fit = glmnet(x, g4, family = "multinomial", alpha = 0, standardize = FALSE, intercept = FALSE, lambda = 0.5)},
#                                IBMR_Newton = {test = fit_Gamma_Newton(list(Y_mat), list(x), list(x), 0.5, 1000, 1e-6, list(matrix(0, ncol = 4, nrow = 100)))},
#                                IBMR = {test = fit_Gamma(list(Y_mat), list(x), list(x), 0.5, 1000, 1e-6, list(matrix(0, ncol = 4, nrow = 100)))}, times = 100)
