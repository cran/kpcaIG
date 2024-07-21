#' Kernel Principal Components Analysis
#'
#' Kernel Principal Components Analysis,  a nonlinear version of principal component analysis obtrained through the so-called kernel trick.
#'

#' @param data The data matrix organized by rows. Users should scale the data appropriately before applying this function, if relevant.
#' @param kernel The kernel function used for the analysis. It can be chosen from the following strings: 'rbfdot': Radial Basis kernel function "Gaussian", 'polydot': Polynomial kernel function, 'vanilladot': Linear kernel function, 'tanhdot': Hyperbolic tangent kernel function.

#' @param kpar The list of hyper-parameters (kernel parameters) used with the kernel function. The valid parameters for each kernel type are as follows:   sigma: inverse kernel width for the Radial Basis kernel function 'rbfdot'. degree, scale, offset for the Polynomial kernel function 'polydot'. scale, offset for the Hyperbolic tangent kernel function 'tanhdot'.
 
#' @param features The number of features (kernel principal components) to use for the analysis. Default: 0 , (all)

#' @return The result of kpca, a Formal class kpca as in library(kernlab). An S4 object containing the principal component vectors along with the corresponding eigenvalues.
#' @importFrom kernlab kpca

#' @export

kernelpca <- function(data, kernel = "vanilladot", kpar = list(), features = 0) {
  requireNamespace("kernlab", quietly = TRUE)
  # KPCA
  kpca_result <- kernlab::kpca(data, kernel = kernel, kpar = kpar, features = features)

  #  The result of kpca
  return(kpca_result)
}
