#' KPCA-IG: variables interpretability in kernel PCA
#'
#' KPCA-IG, kernel pca interpretable gradient. It is the function that gives the feature ranking, from the most to the least relevant variable. The ranking is obtained through the kernel’s partial derivatives computation. A score, which corresponds to the score mean among the sample points, is assigned to each input feature.
#'
#' @param kpca_result The result of the previously obtained kernel PCA analysis.
#' @param dim A vector of kernel principal components to use for the computation of the scores. Each value should be less than or equal to the number of components of the kPCA.
#' @param mean_type Type of mean. Possible values are "arithmetic", "geometric", "harmonic", "trimmed", "median". Default = "arithmetic"
#' @param trim_ratio For mean_type == "trimmed", it is the fraction (0 to 0.5) of scores to be trimmed from each end before the mean is computed for a more robust to outliers arithmetic mean computation
#' @return A data frame containing the sorted variables and their scores sorted in decreasing order.
#' @importFrom kernlab rotated pcv xmatrix kernelMatrix rbfdot
#' @importFrom stats sd
#' @importFrom progress progress_bar
#' @export
kpca_igrad <- function(kpca_result, dim, mean_type = "arithmetic", trim_ratio = 0.1) {
  requireNamespace("kernlab", quietly = TRUE)
  requireNamespace("progress", quietly = TRUE)
  
  # Ensure dim is a vector
  if (!is.vector(dim)) {
    stop("dim should be a vector of component indices.")
  }

  # Determine the number of components in the kpca_result
  ncomp <- ncol(kernlab::rotated(kpca_result))
  
  # Check if the specified components are valid
  if (any(dim > ncomp)) {
    stop("One or more specified components exceed the number of components in the KPCA result.")
  }
  
  p <- ncol(kpca_result@xmatrix)
  mean_norms <- numeric(p)
  std_norms <- numeric(p)

  # Progress bar initialization
  pb <- progress::progress_bar$new(
    format = "  [:bar] :percent in :elapsed, ETA: :eta",
    total = p,
    width = 60
  )

  kernel_type <- class(kpca_result@kernelf)[1]

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

  system.time({
    for (vk in 1:p) {
      gene_var <- colnames(kernlab::xmatrix(kpca_result))[vk]
      xmatrix_col <- kpca_result@xmatrix[, gene_var]
      if (kernel_type == "vanillakernel") {
        dz_dt <- matrix(t(xmatrix_col), nrow = n, ncol = n, byrow = TRUE) %*% MM %*% Vcuc[, dim, drop = FALSE]
      }  else if (kernel_type == "polykernel") {
        K_poly  <- matrix(t(xmatrix_col), nrow = n, ncol = n, byrow = TRUE)
        derive_poly <- degree * scale * K_poly * inner_products
        dz_dt <- derive_poly %*% MM %*% Vcuc[, dim, drop = FALSE]
      } else if (kernel_type == "rbfkernel") {
        K <- outer(xmatrix_col, xmatrix_col, "-")
        dz_dt <- (-2 * sigma * kernel_Dataset * K) %*% MM %*% Vcuc[, dim, drop = FALSE]
      } else if (kernel_type == "tanhkernel") {
        K_tg  <- matrix(t(xmatrix_col), nrow = n, ncol = n, byrow = TRUE)
        derive_tanh <- scale * K_tg * (1 / cosh(scale * inner_products_tanh + offset))
        dz_dt <- derive_tanh %*% MM %*% Vcuc[, dim, drop = FALSE]
      } else {
        stop("Kernel type not supported.")
      }

      norm_squared <- sqrt(rowSums(dz_dt^2))
      
      if (mean_type == "arithmetic") {
        mean_norms[vk] <- mean(norm_squared)
      } else if (mean_type == "geometric") {
        mean_norms[vk] <- exp(mean(log(norm_squared)))
      } else if (mean_type == "harmonic") {
        mean_norms[vk] <- length(norm_squared) / sum(1 / norm_squared)
      } else if (mean_type == "median") {
        mean_norms[vk] <- median(norm_squared)
      } else if (mean_type == "trimmed") {
        mean_norms[vk] <- mean(norm_squared, trim = trim_ratio)
      } else {
        stop("Invalid mean_type. Choose 'arithmetic', 'geometric', 'harmonic', 'median', or 'trimmed'.")
      }

      std_norms[vk] <- stats::sd(norm_squared)

      pb$tick()
    }
  })

  sorted_mean_norms <- sort(mean_norms, decreasing = TRUE)
  indices <- order(mean_norms, decreasing = TRUE)
  sorted_column_names <- colnames(kpca_result@xmatrix)[indices]

  results <- data.frame(
    column_names = sorted_column_names,
    means_norms = sorted_mean_norms,
    std_norms = std_norms[indices]
  )

  return(results)
}

