% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kernelpca.R
\name{kernelpca}
\alias{kernelpca}
\title{Kernel Principal Components Analysis}
\usage{
kernelpca(data, kernel = "vanilladot", kpar = list(), features = 0)
}
\arguments{
\item{data}{The data matrix organized by rows. Users should scale the data appropriately before applying this function, if relevant.}

\item{kernel}{The kernel function used for the analysis. It can be chosen from the following strings:
 \itemize{
    \item{\code{'rbfdot'}}: Radial Basis kernel function "Gaussian"
    \item{\code{'polydot'}}: Polynomial kernel function
    \item{\code{'vanilladot'}}: Linear kernel function
    \item{\code{'tanhdot'}}: Hyperbolic tangent kernel function
  }
}
\item{kpar}{The list of hyper-parameters (kernel parameters) used with the kernel function. The valid parameters for each kernel type are as follows:
\itemize{
  \item{\code{sigma}}: inverse kernel width for the Radial Basis kernel function \code{"rbfdot"}.
  \item{\code{degree}}, \code{scale}, \code{offset} for the Polynomial kernel function \code{"polydot"}.
  \item{\code{scale}}, \code{offset} for the Hyperbolic tangent kernel function \code{"tanhdot"}.
}
}

\item{features}{The number of features (kernel principal components) to use for the analysis. Default: 0 , (all)}
}
\value{
\code{kernelpca} returns an S4 object of formal class kpca as in library(kernlab) containing the principal component vectors along with the corresponding eigenvalues.

\item{pcv}{pcv	a matrix containing the principal component vectors (column wise)}

\item{eig}{The corresponding eigenvalues}

\item{rotated}{The original data projected (rotated) on the principal components}

\item{xmatrix}{The original data matrix}


}
\description{
Kernel Principal Components Analysis,  a nonlinear version of principal component analysis obtrained through the so-called kernel trick.
}


\references{
Scholkopf B., Smola A. and Muller K.R. (1998) Nonlinear component analysis as a kernel eigenvalue problem. Neural Computation, 10, 1299-1319.

}
\examples{
# Example
library(WallomicsData)
library(kpcaIG)
library(ggplot2)
library(kernlab)
Transcriptomics_Stems_s <- scale(Transcriptomics_Stems)

kpca_tan <-  kernelpca(as.matrix(Transcriptomics_Stems_s),
                          kernel = "tanhdot",
                          kpar = list(scale = 0.0001, offset = 0.01))


ggplot(data = data.frame(rotated(kpca_tan), Genetic_Cluster), 
       aes(x = X1, y = X2, shape = Genetic_Cluster)) +
  geom_point(size = 2, aes(color = Genetic_Cluster)) +
  xlab("1st kernel PC") +
  ylab("2nd kernel PC") +
  labs(color = "Genetic_Cluster", shape = "Genetic_Cluster") +
  theme_minimal()


ggplot(data = data.frame(rotated(kpca_tan), Ecotype), 
       aes(x = X1, y = X2, shape =  Ecotype)) +
  geom_point(size = 2, aes(color =  Ecotype)) +
  xlab("1st kernel PC") +
  ylab("2nd kernel PC") +
  labs(color = " Ecotype", shape = " Ecotype") +
  theme_minimal()


#Equivalently with function plot_kpca2D

plot_kpca2D(kpca_tan, groups = Genetic_Cluster)
plot_kpca2D(kpca_tan, groups = Ecotype)


}
