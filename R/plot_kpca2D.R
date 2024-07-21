utils::globalVariables(c("x1", "y1"))

#' 2D Kernel principal analysis plot with variables representation
#'
#' With this function it is possible to visualize an original variable of interest in the first two principal components or any other specified components. The variable is displayed as an arrow, showing its relevance in the relative position of each sample point in the kernel component space.
#'
#' @param kpca_result The result of the previously obtained kernel PCA analysis.
#' @param target_variable A string indicating the name of the variable of interest to visualize as arrows on the kernel PCA plot.
#' @param groups A vector indicating the grouping of data points, if applicable. Default: NULL
#' @param scale Coefficient to adjust the lengths of the arrows. Default: 100
#' @param components A numeric vector of length 2 specifying the indices of the components to plot. Default: c(1, 2)
#' @param arrow_col Colour of the arrows. Default: '#D3D3D3'
#' @param main_title Graph title. Default: "Kernel principal component analysis"
#'
#' @import ggplot2
#' @importFrom viridis scale_fill_viridis_d
#' @export
plot_kpca2D <- function(kpca_result, target_variable, groups = NULL, scale = 100,components = c(1, 2), arrow_col = "#D3D3D3", main_title = "Kernel principal component analysis") {
  requireNamespace("ggplot2", quietly = TRUE)
  
  # Ensure components is a numeric vector of length 2
  if (length(components) != 2 || !all(components %in% 1:ncol(kernlab::rotated(kpca_result)))) {
    stop("Components must be a numeric vector of length 2, and must be valid component indices.")
  }
  
  scores <- data.frame(kernlab::rotated(kpca_result))
  
  if (is.null(groups)) {
    groups <- rep("No group", nrow(kernlab::rotated(kpca_result)))
  }
  
  Groups <- factor(groups)
  
  # Extract the specified components
  x0 <- scores[, components[1]]
  y0 <- scores[, components[2]]
  dz_dt <- igrad_rot(kpca_result, target_variable)
  dz_dt <- data.frame(dz_dt)
  
  # Extract explained variance
  explained_var <- kernlab::eig(kpca_result)
  explained_var <- explained_var / sum(explained_var) * 100  # Convert to percentage
  xlab_text <- sprintf("kPC%d: %.1f%% expl. var", components[1], explained_var[components[1]])
  ylab_text <- sprintf("kPC%d: %.1f%% expl. var", components[2], explained_var[components[2]])
  
  ggplot(scores, aes(x = x0, y = y0, color = Groups, linetype = target_variable)) +
    geom_point(aes(shape = Groups)) +
    xlab(xlab_text) +
    ylab(ylab_text) +
    ggtitle(main_title) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5), # Improved font settings for title
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 10)
    ) +
    geom_segment(data = data.frame(x0 = x0, y0 = y0, x1 = x0 + dz_dt[, components[1]] * scale, y1 = y0 + dz_dt[, components[2]] * scale),
                 aes(x = x0, y = y0, xend = x1, yend = y1),
                 arrow = arrow(length = unit(0.25, "cm")), color = arrow_col) +
    scale_color_viridis_d() +  
    scale_linetype_manual(values = "solid") +
    labs(linetype = "Variable name")
}