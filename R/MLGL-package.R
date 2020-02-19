#' @importFrom graphics abline legend matplot plot par text
#' @importFrom stats anova approx dist lm p.adjust predict binomial glm model.matrix model.response
#' @import gglasso MASS Matrix fastcluster FactoMineR parallelDist
#' 
#' @title MLGL
#' @docType package
#' @aliases MLGL-package
#' @name MLGL-package
#' @description  
#' This package presents a method combining Hierarchical Clustering and Group-lasso. Usually, a single partition of the covariates is used in the group-lasso.
#' Here, we provides several partition from the hierarchical tree.
#' 
#' A post-treatment method based on statistical test (with FWER and FDR control) for selecting the regularization parameter and the optimal group for this value is provided.
#' This method can be applied for the classical group-lasso and our method.  
#'
#' 
#' @details 
#' The function \link{MLGL} performs the hierarchical clustering and the group-lasso. The post-treatment method can be performed with \link{hierarchicalFWER} and \link{selFWER} functions.
#' The whole process can be run with the \link{fullProcess} function.
#' 
#' 
#' @author Quentin Grimonprez 
#' 
#' Maintainer: Quentin Grimonprez  <quentin.grimonprez@@inria.fr>
#'  
#' @references "MLGL: An R package implementing correlated variable selection by hierarchical clustering and group-Lasso.", Quentin Grimonprez, Samuel Blanck, Alain Celisse, Guillemette Marot (2018). \url{https://hal.inria.fr/hal-01857242} 
#' 
#' @examples 
#' # Simulate gaussian data with block-diagonal variance matrix containing 12 blocks of size 5
#' X <- simuBlockGaussian(50, 12, 5, 0.7)
#' # Generate a response variable
#' y <- X[,c(2,7,12)] %*% c(2, 2, -2) + rnorm(50, 0, 0.5)
#' # Apply MLGL method
#' res <- MLGL(X, y)
#' 
#' @seealso \link{MLGL}, \link{cv.MLGL}, \link{fullProcess}, \link{hierarchicalFWER} 
#' 
#' @keywords package
NULL