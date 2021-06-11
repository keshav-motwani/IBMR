# check against CVXR for fine case, coarse case not convex

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

lambda = 0.01
rho = 0.5

system.time({test = fit_alpha_Beta_Gamma(Y_matrix_list, X_list, X_list, lambda, rho, 1000, 1e-6, rep(0, 4), matrix(0, nrow = 20, ncol = 4), lapply(1:length(X_list), function(x) matrix(0, nrow = 20, ncol = 4)))})
system.time({test_Newton = fit_alpha_Beta_Gamma_Newton(Y_matrix_list, X_list, X_list, lambda, rho, 1000, 1e-6, rep(0, 4), matrix(0, nrow = 20, ncol = 4), lapply(1:length(X_list), function(x) matrix(0, nrow = 20, ncol = 4)))})

all(diff(test$objective[test$objective != 0]) <= 0)
all(diff(test_Newton$objective[test_Newton$objective != 0]) <= 0)

thresh <- 1e-24

beta <- Variable(p, 4)
alpha <- Variable(4)
Gamma_1 = Variable(p, 4)
Gamma_2 = Variable(p, 4)

obj = 0
lin_pred_1 = rep(1, nrow(X_list[[1]])) %*% t(alpha) + X_list[[1]] %*% beta + X_list[[1]] %*% Gamma_1
for (i in 1:nrow(X_list[[1]])) {
  obj = obj + log_sum_exp(lin_pred_1[i, ])
  # obj = obj - log(exp(lin_pred_1)[i, ] %*% t(Y_matrix_list[[1]][i, , drop = FALSE])) + log_sum_exp(lin_pred_1[i, ])
}
lin_pred_2 = rep(1, nrow(X_list[[2]])) %*% t(alpha) + X_list[[2]] %*% beta + X_list[[2]] %*% Gamma_2
for (i in 1:nrow(X_list[[2]])) {
  obj = obj + log_sum_exp(lin_pred_2[i, ])
  # obj = obj - log(exp(lin_pred_2)[i, ] %*% t(Y_matrix_list[[2]][i, , drop = FALSE])) + log_sum_exp(lin_pred_2[i, ])
}
obj = obj - CVXR::matrix_trace(t(Y_matrix_list[[1]]) %*% lin_pred_1) - CVXR::matrix_trace(t(Y_matrix_list[[2]]) %*% lin_pred_2)
obj = obj / (nrow(X_list[[1]]) + nrow(X_list[[2]]))
obj = obj + lambda * sum(norm2(beta, axis = 1)) + (rho * sum(Gamma_1 ^ 2) / 2) + (rho * sum(Gamma_2 ^ 2) / 2)

prob <- Problem(Minimize(obj))
result <- solve(prob, FEASTOL = thresh, RELTOL = thresh, ABSTOL = thresh, verbose = TRUE, num_iter = 1000)

par(mfrow = c(2, 3))

plot(result$getValue(beta), test$Beta)
abline(0, 1)

test$Gamma_list[[1]] = matrix(test$Gamma_list[[1]], ncol = 4)
plot(result$getValue(Gamma_1), test$Gamma_list[[1]])
abline(0, 1)

test$Gamma_list[[2]] = matrix(test$Gamma_list[[2]], ncol = 4)
plot(result$getValue(Gamma_2), test$Gamma_list[[2]])
abline(0, 1)

plot(result$getValue(beta), test_Newton$Beta)
abline(0, 1)

test_Newton$Gamma_list[[1]] = matrix(test_Newton$Gamma_list[[1]], ncol = 4)
plot(result$getValue(Gamma_1), test_Newton$Gamma_list[[1]])
abline(0, 1)

test_Newton$Gamma_list[[2]] = matrix(test_Newton$Gamma_list[[2]], ncol = 4)
plot(result$getValue(Gamma_2), test_Newton$Gamma_list[[2]])
abline(0, 1)

result$getValue(alpha) - mean(result$getValue(alpha))
test_Newton$alpha

rm(list = ls())
