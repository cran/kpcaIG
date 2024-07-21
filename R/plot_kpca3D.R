#' 3D Kernel principal analysis plot with variables representation
#' @param kpca_result The result of the previously obtained kernel PCA analysis.
#' @param target_variable A string indicating the name of the variable to visualize as arrows on the kernel PCA plot.
#' @param groups A vector indicating the grouping of data points, if applicable. Default: NULL
#' @param scale Coefficient to adjust the lengths of the arrows. Default: 1
#' @param type A character indicating the type of point for the observations. Supported types are: 'p' for points, 's' for spheres. Default: 's'
#' @param size The size of the plotted points. Default: 3/4
#' @param arrow_col Colour of the arrows. Default: '#999999'
#' @param angles Number of barbs of the arrows. Default: 12
#' @param main Graph title. Default: NULL
#' @return Returns an interactive 3D plot.
#' @import rgl
#' @import grDevices
#' @import viridis

#' @export
plot_kpca3D <- function(kpca_result, target_variable, groups = NULL, scale=1, type = "s", size = 3/4, arrow_col = "#999999", angles = 12, main = NULL ) {
  requireNamespace("rgl", quietly = TRUE)
  requireNamespace("grDevices", quietly = TRUE)
  requireNamespace("kernlab", quietly = TRUE)
  requireNamespace("viridis", quietly = TRUE)
  
  if (ncol(kernlab::rotated(kpca_result)) < 3) {
    stop("The kpca_result object must have at least 3 principal components for 3D plotting.")
  }
  
  
  
  dz_dt <- igrad_rot(kpca_result,target_variable)


  # Default to a single color if col is not provided
  if (is.null(groups)) {
    groups <- rep(1, nrow(kernlab::rotated(kpca_result)))  # Default to one color for all points
  }
  
colors <- viridis::viridis(length(unique(groups)), option = "D")
groups <- colors[as.numeric(factor(groups))]

  # Extraction of kenrel principal components
  x0 <- kernlab::rotated(kpca_result)[, 1]
  y0 <- kernlab::rotated(kpca_result)[, 2]
  z0 <- kernlab::rotated(kpca_result)[, 3] 

  # rgl window
    rgl::open3d()
    rgl::clear3d()
    rgl::bg3d("#fbfeff")

  # Tip of the arrows
  x1 <- x0 + scale * dz_dt[, 1]
  y1 <- y0 + scale * dz_dt[, 2]
  z1 <- z0 + scale * dz_dt[, 3] 
 
  
    rgl::plot3d(x = x0, y = y0, z = z0,
                xlab = "", ylab = "", zlab = "",
                col = groups, type = type, size = size, decorate = TRUE) 
 
 rgl::decorate3d(main = main,expand = 1.2, box =  FALSE, xlab = "kPCA1",  ylab = "kPCA2",  zlab = "kPCA3")
 
  n_points <- length(x0)

  for (i in 1:n_points) {
    p0 <- c(x0[i], y0[i], z0[i])
    p1 <- c(x1[i], y1[i], z1[i])
    
   
    rgl::arrow3d(p0 = p0, p1 = p1, type = "rotation", col = arrow_col, n = angles,
                 width = 0.3, s = 0.2 , theta = 0.2 )
  }

  # Graph
    return(rgl::rglwidget())
  
}
