#' Internal function for the computation of the gradient and coordinates of data points into the kpca axes for a SINGLE variable 
#'
#' @param kpca_result Result from the kernel PCA analysis.
#' @param target_variable The variable of interest for gradient calculation.
#' @return A data frame containing the derivatives with respect to the variable of interest projected into the kernel axes
#' @keywords internal

igrad_rot <- function(kpca_result, target_variable) {
  requireNamespace("kernlab", quietly = TRUE)
  requireNamespace("progress", quietly = TRUE)
  kernel_type <- class(kpca_result@kernelf)[1]
  xmatrix_col <- kpca_result@xmatrix[, target_variable]
  sigma <- if (kernel_type == "rbfkernel") kpca_result@kernelf@kpar$sigma else NULL
  degree <- if (kernel_type == "polykernel") kpca_result@kernelf@kpar$degree else NULL
  scale <- if (kernel_type %in% c("polykernel", "tanhkernel")) kpca_result@kernelf@kpar$scale else NULL
  offset <- if (kernel_type %in% c("polykernel", "tanhkernel")) kpca_result@kernelf@kpar$offset else NULL



  if (kernel_type == "polykernel") {
    inner_products <- (scale * (kpca_result@xmatrix %*% t(kpca_result@xmatrix)) + offset)^(degree - 1)
  } else if (kernel_type == "tanhkernel") {
    inner_products_tanh <- kpca_result@xmatrix %*% t(kpca_result@xmatrix)
  }

  kernel_Dataset <- if (kernel_type == "rbfkernel") kernlab::kernelMatrix(kernlab::rbfdot(sigma), as.matrix(kpca_result@xmatrix)) else NULL

  n <- nrow(kpca_result@xmatrix)
  MM <- diag(n) - matrix(1/n, n, n)
  Vcuc <- kernlab::pcv(kpca_result)



  if (kernel_type == "vanillakernel") {
    dz_dt <- matrix(t(xmatrix_col), nrow = n, ncol = n, byrow = TRUE) %*% MM %*% Vcuc
  }  else if (kernel_type == "polykernel") {
    K_poly=matrix(t(xmatrix_col), nrow = n, ncol = n, byrow = TRUE)
    derive_poly <- degree * scale * K_poly * inner_products
    dz_dt <- derive_poly %*% MM %*% Vcuc
  } else if (kernel_type == "rbfkernel") {
    K <- outer(xmatrix_col, xmatrix_col, "-")
    dz_dt <- (-2 * sigma * kernel_Dataset * K) %*% MM %*% Vcuc
  } else if (kernel_type == "tanhkernel") {
    K_tg=matrix(t(xmatrix_col), nrow = n, ncol = n, byrow = TRUE)
    derive_tanh <- scale*K_tg * (1/cosh(scale * inner_products_tanh + offset))
    dz_dt <- derive_tanh %*% MM %*% Vcuc
  } else {
    stop("Kernel type not supported.")
  }

  return(dz_dt)
}
