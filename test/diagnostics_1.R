# check that fitted Gamma of fit_Gamma is same as ridge regression when only 1 dataset

set.seed(1)

library(IBMR)
library(glmnet)

create_Y_matrix = function(Y, categories) {

  Y = t(sapply(Y, function(y) categories == y) * 1)

  colnames(Y) = categories

  return(Y)

}

n = 10000

g4 = sample(1:4, n, replace = TRUE)
x = matrix(rnorm(n * 100), n, 100)

system.time({fit3 = glmnet(x, g4, family = "multinomial", alpha = 0, standardize = FALSE, intercept = FALSE, lambda = 0.5)})
test2 = as.matrix(do.call(cbind, coef(fit3, 0))[-1, ])

Y_mat = create_Y_matrix(g4, 1:4)

system.time({test = fit_Gamma(list(Y_mat), list(x), list(x), 0.5, 1000, 1e-6, list(matrix(0, ncol = 4, nrow = 100)))$Gamma_list})
test = matrix(test[[1]], ncol = 4)
plot(test[, 1] - test[, 4], test2[, 1] - test2[, 4])
abline(0, 1)

microbenchmark::microbenchmark(glmnet = {fit3 = glmnet(x, g4, family = "multinomial", alpha = 0, standardize = FALSE, intercept = FALSE, lambda = 0.5)},
                               IBMR_Newton = {test = fit_Gamma_Newton(list(Y_mat), list(x), list(x), 0.5, 1000, 1e-6, list(matrix(0, ncol = 4, nrow = 100)))},
                               IBMR = {test = fit_Gamma(list(Y_mat), list(x), list(x), 0.5, 1000, 1e-6, list(matrix(0, ncol = 4, nrow = 100)))}, times = 100)

rm(list = ls())
