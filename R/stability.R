#' Stability selection for \code{\link{MLGL}}
#'
#' @title Stability Selection for Multi-Layer Group-lasso
#'
#' @author Quentin Grimonprez
#' @param X matrix of size n*p
#' @param y vector of size n. If loss = "logit", elements of y must be in {-1,1} 
#' @param B number of bootstrap sample
#' @param fraction Fraction of data used at each of the \code{B} sub-samples
#' @param lambda lambda values for group lasso. If not provided, the function generates its own values of lambda
#' @param hc output of \code{\link{hclust}} function. If not provided, \code{\link{hclust}} is run with ward.D2 method
#' @param weightLevel a vector of size p for each level of the hierarchy. A zero indicates that the level will be ignored. If not provided, use 1/(height between 2 successive levels)
#' @param weightSizeGroup a vector
#' @param loss a character string specifying the loss function to use, valid options are: "ls" least squares loss (regression) and "logit" logistic loss (classification)
#' @param intercept should an intercept be included in the model ?
#' @param verbose print some informations
#' @param ... Others parameters for \code{\link{gglasso}} function
#'
#' @return a stability.MLGL object containing :
#' \describe{
#' \item{lambda}{sequence of \code{lambda}.}
#' \item{B}{Number of bootstrap samples.}
#' \item{stability}{A matrix of size length(lambda)*number of groups containing the probability of selection of each group}
#' \item{var}{vector containing the index of covariates}
#' \item{group}{vector containing the index of associated groups of covariates}
#' \item{time}{computation time}
#' } 
#'
#' @references Meinshausen and Buhlmann (2010). Stability selection. In : Journal of the Royal Statistical Society : Series B (Statistical Methodology) 72.4, p. 417-473.
#'
#' @details 
#' Hierarhical clustering is performed with all the variables. Then, the partitions from the different
#' levels of the hierarchy are used in the differents run of MLGL for estimating the probability of selection of each group.
#'  
#' @examples 
#' \donttest{
#' set.seed(42)
#' # Simulate gaussian data with block-diagonal variance matrix containing 12 blocks of size 5
#' X <- simuBlockGaussian(50, 12, 5, 0.7)
#' 
#' # Generate a response variable
#' y <- X[,c(2,7,12)]%*%c(2,2,-2) + rnorm(50, 0, 0.5)
#' 
#' # Apply stability.MLGL method
#' res <- stability.MLGL(X, y)
#' }
#' 
#' @seealso \link{cv.MLGL}, \link{MLGL}
#' 
#' @export
stability.MLGL <- function(X, y, B = 50, fraction = 0.5, hc = NULL, lambda = NULL, weightLevel = NULL, weightSizeGroup = NULL, loss = c("ls", "logit"), intercept = TRUE, verbose = FALSE,...)
{
  
  #check parameters 
  loss = match.arg(loss)
  .checkParameters(X, y, hc, lambda, weightLevel, weightSizeGroup, intercept, verbose, loss, NULL)
  
  #B
  if( !is.numeric(B) | (length(B)!=1) )
    stop("B must be a positive integer.")
  if(!.is.wholenumber(B)) 
    stop("B must be a positive integer.")
  #fraction
  if( !is.numeric(fraction) | (length(fraction)!=1) )
    stop("fraction must be a positive real lesser than 1.")
  if( (fraction<0) | (fraction>1) ) 
    stop("fraction must be a positive real lesser than 1.")
  
  ################################ same as MLGL function
  #define some usefull variables
  n <- nrow(X)
  p <- ncol(X)
  tcah <- rep(NA,3)
  
  #if no hc output provided, we make one
  if(is.null(hc))
  {
    if(verbose)
      cat("Computing hierarchical clustering...")
    
    t1 <- proc.time()
    d <- dist(t(X))
    hc <- hclust(d,method="ward.D2")
    t2 <- proc.time()
    tcah <- t2-t1
    if(verbose)
      cat("DONE in ",tcah[3],"s\n")
  }
  
  #compute weight, active variables and groups
  if(verbose)
    cat("Preliminary step...")
  t1 = proc.time()
  prelim <- preliminaryStep(hc,weightLevel,weightSizeGroup)  
  
  #duplicate data
  Xb <- X[,prelim$var]
  t2 = proc.time()
  if(verbose)
    cat("DONE in ",(t2-t1)[3],"s\n")
  
  ################################ END same as MLGL function
  
  
  ######## stability selection
  if(verbose)
    cat("Running stability selection...\n")
  t1 = proc.time()
  
  #create lambda sequence, if not provided
  if(is.null(lambda))
    lambda = lambdaseq(Xb, y, prelim$group, prelim$weight, intercept = intercept)
  
  proba <- Matrix(0, nrow = length(lambda), ncol = prelim$group[length(prelim$group)])  
  
  #bootstrap
  for(b in 1:B)
  {
    t1b <- proc.time()
    if(verbose)
      cat("   Running sample",b,"...")
    #sample of size n/2
    testind <- sample(1:n,floor(n*fraction))
    
    t2 = proc.time()
    res <- gglasso(Xb[testind,], y[testind], prelim$group, lambda = lambda, pf = prelim$weight, intercept = intercept, loss = loss, ...)
    t3 <- proc.time()
    
    #record selected groups
    for(indlam in 1:length(lambda))
    {
      non0 <- which(res$beta[,indlam]!=0)#non 0 coefficients
      non0gr <- unique(prelim$group[non0])#associated non 0 groups
      proba[indlam, non0gr] = proba[indlam,non0gr]+1 #increment counter of selected groups
    }#fin for lambda
    t2b <- proc.time()
    
    if(verbose)
      cat("DONE in",(t2b-t1b)[3],"s \r")
  }#fin for B
  
  #divide by B for having probabilities
  proba = proba/B
  t2 = proc.time()
  
  if(verbose)
    cat("\nStability selection DONE in",(t2-t1)[3],"s")
  
  #create output
  res <- list()
  res$lambda = lambda
  res$stability = proba
  res$B = B
  res$var = prelim$var
  res$group = prelim$group
  res$time = c(tcah[3],(t2-t1)[3])
  names(res$time) = c("hclust","stability")   
  class(res) = "stability.MLGL"
  
  
  return(res)
}


#
# generate lambda sequence for group lasso
#
# @param X data matrix
# @param y response
# @param group vector indicates the grouping of covariates
# @param lambda.factor the ratio lambdamin/lambdamax
# @param length length of the lambda sequence
#
# @return lambda sequence
#
lambdaseq <- function(X, y, group, weight, lambda.factor = NULL, length = 100, intercept = TRUE)
{
  #dimension of X
  n <- nrow(X)
  p <- ncol(X)
  
  #compute the max lambda by KKT
  if(intercept)
    y = y-mean(y)
  l <- tapply(1:p,group,FUN=function(x){sqrt(sum((t(X[,x])%*%y)^2))/n})
  l = l/weight
  lambdamax <- max(l)
  
  #heuristic min lambda
  if(is.null(lambda.factor))
    lambda.factor = ifelse( n < p, 0.05, 0.001)
  
  lambdamin <- lambdamax*lambda.factor
  
  #generate sequence on a log scale
  lambda <- exp(seq(log(lambdamax),log(lambdamin),length=length))

  return(lambda)  
}

