#' @export
check_KKT_IBMR = function(Y_matrix_list,
                          X_list,
                          Z_list,
                          lambda,
                          rho,
                          alpha,
                          Beta,
                          Gamma_list) {

  N = sum(sapply(Y_matrix_list, nrow))

  grad_alpha = IBMR:::compute_gradient_alpha(Y_matrix_list, X_list, Z_list, alpha, Beta, Gamma_list, N)
  KKT_alpha = max(abs(grad_alpha))

  nonzero = rowSums(Beta) != 0
  grad_Beta = IBMR:::compute_gradient_Beta(Y_matrix_list, X_list, Z_list, alpha, Beta, Gamma_list, N)
  KKT_Beta = suppressWarnings(max(abs(grad_Beta[nonzero, , drop = FALSE] + lambda * (Beta[nonzero, , drop = FALSE] / apply(Beta[nonzero, , drop = FALSE], 1, function(x) sqrt(sum(x ^ 2)))))))
  KKT_Beta_0 = suppressWarnings(max(apply(-grad_Beta[!nonzero, , drop = FALSE] / lambda, 1, function(x) sqrt(sum(x ^ 2)))))
  if (sum(nonzero) == 0) KKT_Beta = 0
  if (sum(!nonzero) == 0) KKT_Beta_0 = 1

  KKT_Gamma = max(sapply(1:length(Y_matrix_list), function(i) abs(max(IBMR:::compute_gradient_Gamma(Y_matrix_list[[i]], X_list[[i]], Z_list[[i]], alpha, Beta, Gamma_list[[i]], rho, N)))))

  return(c(KKT_alpha = KKT_alpha, KKT_Beta = KKT_Beta, KKT_Beta_0 = KKT_Beta_0, KKT_Gamma = KKT_Gamma))

}

#' @export
check_KKT_IBMR_no_Gamma = function(Y_matrix_list,
                                   X_list,
                                   lambda,
                                   alpha,
                                   Beta) {

  Z_list = lapply(Y_matrix_list, function(Y) matrix(1, nrow = nrow(Y), ncol = 1))
  Gamma_list = lapply(Z_list, function(x) matrix(0, nrow = ncol(x), ncol = ncol(Y_matrix_list[[1]])))

  N = sum(sapply(Y_matrix_list, nrow))

  grad_alpha = IBMR:::compute_gradient_alpha(Y_matrix_list, X_list, Z_list, alpha, Beta, Gamma_list, N)
  KKT_alpha = max(abs(grad_alpha))

  nonzero = rowSums(Beta) != 0
  grad_Beta = IBMR:::compute_gradient_Beta(Y_matrix_list, X_list, Z_list, alpha, Beta, Gamma_list, N)
  KKT_Beta = suppressWarnings(max(abs(grad_Beta[nonzero, , drop = FALSE] + lambda * (Beta[nonzero, , drop = FALSE] / apply(Beta[nonzero, , drop = FALSE], 1, function(x) sqrt(sum(x ^ 2)))))))
  KKT_Beta_0 = suppressWarnings(max(apply(-grad_Beta[!nonzero, , drop = FALSE] / lambda, 1, function(x) sqrt(sum(x ^ 2)))))
  if (sum(nonzero) == 0) KKT_Beta = 0
  if (sum(!nonzero) == 0) KKT_Beta_0 = 1

  return(c(KKT_alpha = KKT_alpha, KKT_Beta = KKT_Beta, KKT_Beta_0 = KKT_Beta_0))

}
