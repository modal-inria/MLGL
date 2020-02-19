#'
#' Group-lasso with overlapping groups
#' 
#' @param X matrix of size n*p
#' @param y vector of size n. If loss = "logit", elements of y must be in {-1,1} 
#' @param var vector containing the variable to use
#' @param group vector containing the associated groups
#' @param lambda lambda values for group lasso. If not provided, the function generates its own values of lambda
#' @param weight a vector the weight for each group. Default is the sqaure root of the size of each group
#' @param loss a character string specifying the loss function to use, valid options are: "ls" least squares loss (regression) and "logit" logistic loss (classification)
#' @param intercept should an intercept be included in the model ?
#' @param ... Others parameters for \code{\link{gglasso}} function
#' 
#' @return a MLGL object containing :
#' \describe{
#'   \item{lambda}{lambda values}
#'   \item{b0}{intercept values for \code{lambda}}
#'   \item{beta}{A list containing the values of estimated coefficients for each values of \code{lambda}}
#'   \item{var}{A list containing the index of selected variables for each values of \code{lambda}}
#'   \item{group}{A list containing the values index of selected groups for each values of \code{lambda}}
#'   \item{nVar}{A vector containing the number of non zero coefficients for each values of \code{lambda}}
#'   \item{nGroup}{A vector containing the number of non zero groups for each values of \code{lambda}}
#'   \item{structure}{A list containing 3 vectors. var : all variables used. group : associated groups. 
#'   weight : weight associated with the different groups.}
#'   \item{time}{computation time}
#'   \item{dim}{dimension of \code{X}}
#' } 
#'
#' @details Use a group-lasso algorithm (see \code{\link{gglasso}}) to solve a group-lasso with overlapping groups. 
#' Each variable j of the original matrix \code{X} is paste k(j) times in a new dataset with k(j) the number of different groups containing the variable j.
#' The new dataset is used to solve the group-lasso with overlapping groups running a group-lasso algorithm. 
#'
#'
#' @source Laurent Jacob, Guillaume Obozinski, and Jean-Philippe Vert. 2009. Group lasso with overlap and graph lasso. In Proceedings of the 26th Annual International Conference on Machine Learning (ICML '09).
#'
#'
#' @examples 
#' # Least square loss
#' set.seed(42)
#' X <- simuBlockGaussian(50, 12, 5, 0.7)
#' y <- X[, c(2, 7, 12)]%*%c(2, 2, -2) + rnorm(50, 0, 0.5)
#' var <- c(1:60, 1:8, 7:15)
#' group <- c(rep(1:12, each = 5), rep(13, 8), rep(14, 9))
#' res <- overlapgglasso(X, y, var, group)
#'
#'# Logistic loss
#' y <- 2*(rowSums(X[,1:4])>0)-1
#' var <- c(1:60, 1:8, 7:15)
#' group <- c(rep(1:12, each = 5), rep(13, 8), rep(14, 9))
#' res <- overlapgglasso(X, y, var, group, loss = "logit")
#'
#'
#'
#'
#' @seealso \code{\link{listToMatrix}}
#'
#' @export
overlapgglasso <- function(X, y, var, group, lambda = NULL, weight = NULL, loss = c("ls", "logit"), intercept = TRUE,...)
{
  #check parameters
  .checkOverlap(X, y, var, group, lambda, weight, intercept)
  loss <- match.arg(loss)
  
  #order group (for gglasso)
  ord <- order(group)
  groupord <- group[ord]
  #order var according to group
  varord <- var[ord]
  
  #transform group to have consecutive numbers (for gglasso)
  groupb <- cumsum(!duplicated(groupord))
  
  #new data
  Xb <- X[,varord]
  
  #if weight not provided, use the sqrt of the size of groups
  if(is.null(weight))
    weight = as.numeric(sqrt(table(groupb)))
  
  #overlap group lasso
  t1 <- proc.time()
  res <- gglasso(Xb, y, groupb, pf = weight, lambda = lambda, intercept = intercept, loss = loss,...)
  t2 <- proc.time()
  

  #create output object
  res2 <- list()
  res2$lambda = res$lambda
  non0 = apply(res$beta,2, FUN=function(x){which(x!=0)})
  res2$var = lapply(non0, FUN=function(x){varord[x]})
  res2$nVar = sapply(res2$var, FUN=function(x){length(unique(x))})
  res2$group = lapply(non0, FUN=function(x){groupb[x]})
  res2$nGroup = sapply(res2$group, FUN=function(x){length(unique(x))})
  res2$beta = lapply(1:length(res$lambda), FUN=function(x){res$beta[non0[[x]],x]})
  res2$b0 = res$b0
  res2$structure = list(group=groupb, var=varord, weight=weight)
  res2$dim = dim(X)
  res2$time = (t2-t1)[3]    
  class(res2) = "MLGL"
  
  return(res2)
}




#check parameters for overlapgglasso function
.checkOverlap <- function(X, y, var, group, lambda = NULL, weight = NULL, intercept = TRUE)
{
  #check X
  if(!is.matrix(X)) 
    stop("X has to be a matrix.")
  if(any(is.na(X))) 
    stop("Missing values in X not allowed.")
  if(!is.numeric(X))
    stop("X has to be a matrix of real.")
  
  #check y
  if(!is.numeric(y))
    stop("y has to be a matrix of real.")
  if(any(is.na(y))) 
    stop("Missing values in y not allowed.")
  
  #check if X and y are compatible
  if(nrow(X)!=length(drop(y)))
    stop("The length of y and the number of rows of X don't match.")
  
  #check if var is a vector of positive integer
  if(!is.null(var))
  {
    if(!is.numeric(var))
      stop("var must be a vector of positive integer.")
    if(!all(.is.wholenumber(var)))
      stop("var must be a vector of positive integer.")
    if(any(var < 0))
      stop("var must be a vector of positive integer.")
    if(any(var > ncol(X)))
      stop("An element in var is greater that the number of variables")
  }
  
  #check if group is a vector of positive integer
  if(!is.null(group))
  {
    if(!is.numeric(group))
      stop("group must be a vector of positive integer.")
    if(!all(.is.wholenumber(group)))
      stop("group must be a vector of positive integer.")
    if(any(group < 0))
      stop("group must be a vector of positive integer.")
    if(length(var)!=length(group))
      stop("var and group must have the same size.")
  }
  
  #check if lambda is a vector of positive real
  if(!is.null(lambda))
  {
    if(!is.numeric(lambda))
      stop("lambda must be a vector of positive real.")
    if(any(lambda < 0))
      stop("lambda must be a vector of positive real.")
  }
  
  #check if weight is a vector of positive real
  if(!is.null(weight))
  {
    if(!is.numeric(weight))
      stop("weight must be a vector of positive real.")
    if(length(weight)!=sum(!duplicated(weight)))
      stop("weight must have the same length as the number of differents groups.")
    if(any(weight < 0))
      stop("weight must be a vector of positive real.")
  }
  
  #check if intercept is a boolean
  if(length(intercept)!=1)
    stop("intercept must be a boolean.")
  if(!is.logical(intercept))
    stop("intercept must be a boolean.")
  
  invisible(return(NULL))
}