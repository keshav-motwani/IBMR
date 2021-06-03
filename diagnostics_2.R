# check that fitted Gamma is same as group lasso regression when only 1 dataset

set.seed(1)

library(IBMR)
library(glmnet)

create_Y_matrix = function(Y, categories) {

  Y = t(sapply(Y, function(y) categories == y) * 1)

  colnames(Y) = categories

  return(Y)

}

n = 1000

g4 = sample(1:4, n, replace = TRUE)
x = matrix(rnorm(n * 100), n, 100)

system.time({fit3 = glmnet(x, g4, family = "multinomial", alpha = 1, standardize = FALSE, type.multinomial = "grouped", maxit = 1e7, thresh = 1e-12)})
test2 = as.matrix(do.call(cbind, coef(fit3, fit3$lambda[20]))[-1, ])

Y_mat = create_Y_matrix(g4, 1:4)

system.time({test = fit_alpha_Beta(list(Y_mat), list(x), list(x), fit3$lambda[20], 1000, 1e-6, rep(0, 4), matrix(0, nrow = 100, ncol = 4))$Beta})
plot(test[, 1] - test[, 4], test2[, 1] - test2[, 4])
abline(0, 1)
