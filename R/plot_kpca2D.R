utils::globalVariables(c("x1", "y1"))

#' 2D Kernel principal analysis plot with variables representation
#'
#' With this function it is possible to visualize an original variable of interest in the first two principal components or any other specified components. The variable is displayed as an arrow, showing its relevance in the relative position of each sample point in the kernel component space.
#'
#' @param kpca_result The result of the previously obtained kernel PCA analysis.
#' @param target_variable A string indicating the name of the variable of interest to visualize as arrows on the kernel PCA plot. Default: NULL
#' @param groups A vector indicating the grouping of data points, if applicable. Default: NULL
#' @param scale Coefficient to adjust the lengths of the arrows. Default: 100
#' @param components A numeric vector of length 2 specifying the indices of the components to plot. Default: c(1, 2)
#' @param arrow_col Colour of the arrows. Default: '#D3D3D3'
#' @param main_title Graph title. Default: "Kernel principal component analysis"
#' @param point_size Sample points size. Default: 2
#' @param arrow_thickness Arrows size. Default: 0.5
#' @param text_size  Text size. Default: 16
#' @param legend_text_size Legend text size. Default: 11
#' @param axis_text_size Axes text size. Default: 12

#'
#' @import ggplot2
#' @importFrom utils modifyList
#' @importFrom viridis scale_fill_viridis_d
#' @export


plot_kpca2D <- function(kpca_result, target_variable = NULL, groups = NULL, scale = 100, components = c(1, 2), arrow_col = "#D3D3D3", main_title = "Kernel principal component analysis", point_size = 2, arrow_thickness = 0.5, text_size = 16, legend_text_size = 11, axis_text_size = 12) {
  requireNamespace("ggplot2", quietly = TRUE)
  
  # Ensure components is a numeric vector of length 2 and valid indices
  if (length(components) != 2 || !all(components %in% 1:ncol(kernlab::rotated(kpca_result)))) {
    stop("Components must be a numeric vector of length 2, and must be valid component indices.")
  }
  
  scores <- data.frame(kernlab::rotated(kpca_result))
  
  if (is.null(groups)) {
    groups <- rep("No group", nrow(scores))
  }
  Groups <- factor(groups)
  
  # Extract the specified components from the scores
  x0 <- scores[, components[1]]
  y0 <- scores[, components[2]]
  
  # Extract explained variance and create axis labels
  explained_var <- kernlab::eig(kpca_result)
  explained_var <- explained_var / sum(explained_var) * 100  # Convert to percentage
  xlab_text <- sprintf("kPC%d: %.1f%% expl. var", components[1], explained_var[components[1]])
  ylab_text <- sprintf("kPC%d: %.1f%% expl. var", components[2], explained_var[components[2]])
  
  # Set up the base ggplot mapping
  base_mapping <- aes(x = x0, y = y0, color = Groups, shape = Groups)
  
  # If target_variable is provided, include a linetype mapping; otherwise, do not.
  if (!is.null(target_variable)) {
    base_mapping <- modifyList(base_mapping, aes(linetype = target_variable))
  }
  
  p <- ggplot(scores, base_mapping) +
    geom_point(size = point_size) +
    xlab(xlab_text) +
    ylab(ylab_text) +
    ggtitle(main_title) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = text_size, face = "bold", hjust = 0.5),
      legend.title = element_text(size = legend_text_size, face = "bold"),
      legend.text = element_text(size = legend_text_size),
      axis.title.x = element_text(size = axis_text_size),
      axis.title.y = element_text(size = axis_text_size)
    )
  
  # If target_variable is provided, compute gradient and add arrows
  if (!is.null(target_variable)) {
    dz_dt <- igrad_rot(kpca_result, target_variable)
    dz_dt <- data.frame(dz_dt)
    
    arrow_data <- data.frame(
      x0 = x0,
      y0 = y0,
      x1 = x0 + dz_dt[, components[1]] * scale,
      y1 = y0 + dz_dt[, components[2]] * scale
    )
    
    p <- p +
      geom_segment(
        data = arrow_data,
        aes(x = x0, y = y0, xend = x1, yend = y1),
        arrow = arrow(length = unit(0.25, "cm")),
        color = arrow_col,
        size = arrow_thickness
      ) +
      scale_linetype_manual(values = "solid") +
      labs(linetype = "Variable name")
  }
  
  p + scale_color_viridis_d()
}