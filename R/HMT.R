#' @title Hierarchical Multiple Testing procedure
#'
#' Apply Hierarchical Multiple Testing procedure on a \code{\link{MLGL}} object
#'
#'  
#' @param res \code{\link{MLGL}} object
#' @param X matrix of size n*p
#' @param y vector of size n.
#' @param control either "FDR" or "FWER"
#' @param alpha control level for testing procedure
#' @param test test used in the testing procedure. Default is \link{partialFtest}
#' @param ... extra parameters fpr \link{selFDR}
#' 
#' @return a list containing :
#' \describe{
#'   \item{lambdaOpt}{lambda values maximizing the number of rejects}
#'   \item{var}{A vector containing the index of selected variables for the first \code{lambdaOpt} value}
#'   \item{group}{A vector containing the values index of selected groups for the first \code{lambdaOpt} value}
#'   \item{selectedGroups}{Selected groups for the first \code{lambdaOpt} value}
#'   \item{reject}{Selected groups for all lambda values}
#'   \item{alpha}{Control level}
#'   \item{test}{Test used in the testing procedure}
#'   \item{control}{"FDR" or "FWER"}
#'   \item{time}{Elapsed time}
#'   \item{hierTest}{list containing the output of the testing function for each lambda. Each element can be used with the \link{selFWER} or \link{selFDR} functions.}
#'   \item{lambda}{lambda path}
#'   \item{nGroup}{Number of groups before testing}
#'   \item{nSelectedGroup}{Numer of groups after testing}
#' } 
#' 
#' @examples
#' set.seed(42)
#' X <- simuBlockGaussian(50, 12, 5, 0.7)
#' y <- X[,c(2,7,12)] %*% c(2,2,-2) + rnorm(50, 0, 0.5)
#' res <- MLGL(X, y)
#' 
#' # perform hierarchical testing with FWER control
#' out <- HMT(res, X, y, alpha = 0.05)
#' 
#' # test a new value of alpha for a specific lambda
#' selFWER(out$hierTest[[60]], alpha = 0.1)
#' 
#' 
#' @seealso \link{hierarchicalFWER} \link{hierarchicalFDR} \link{selFWER} \link{selFDR}
#' 
#' @export
HMT <- function(res, X, y, control = c("FWER", "FDR"), alpha = 0.05, test = partialFtest, ...)
{
  control = match.arg(control)
  
  t1 <- proc.time()
  # choose the right function
  hierTestFunction <- hierarchicalFWER
  selFunction <- selFWER
  if(control == "FDR")
  {
    hierTestFunction = hierarchicalFDR
    selFunction = selFDR  
  }
  
  # testing procedure for each lambda
  REJECT = TEST <- list()
  nbReject <- rep(0, length(res$lambda))
  prevSelGroup = selGroup <- c()
  for(i in 1:length(res$lambda))
  {
    # if no groups are selected we do nothing
    if(length(res$group[[i]])>0)
    {
      
      selGroup = unique(res$group[[i]])
      
      # if the selected groups have not changed compared with the last iteration, we copy the result
      if(setequal(prevSelGroup, selGroup))
      {
        TEST[[i]] = TEST[[i-1]]
        REJECT[[i]] = REJECT[[i-1]]
        nbReject[i] = nbReject[i-1]
      }
      else
      {
        # hierarchical testing and selection
        resTest <- hierTestFunction(X, y, res$group[[i]], res$var[[i]], test)
        
        resSel <- selFunction(resTest, alpha, ...)
        
        # keep outerNode (need for FDR outer = FALSE, do not change in other cases)
        groupSel <- outerNode(resSel$toSel, resTest$hierMatrix)
        # Id of rejected groups
        REJECT[[i]] = (resSel$groupId[resSel$toSel])[groupSel] 
        TEST[[i]] = resTest
        
        # number of rejects for the lambda value
        nbReject[i] = length(REJECT[[i]])
      }
      
    }# end if no selection
    
    prevSelGroup = selGroup
    
  }# end for lambda
  t2 <- proc.time()
  
  hierTestTime <- (t2-t1)[3]
  time = c(res$time, hierTestTime)
  names(time) = c(names(res$time), "test")
  
  # indice of optimal lambda : the one with the greatest number of reject
  indLambdaOpt <- which(nbReject == max(nbReject))
  
  
  
  # group selected for the lambda optimal
  indGroupSel <- res$group[[indLambdaOpt[1]]] %in% REJECT[[indLambdaOpt[1]]]
  group <- res$group[[indLambdaOpt[1]]][indGroupSel]
  var   <- res$var[[indLambdaOpt[1]]][indGroupSel]
  
  
  out <- list(lambdaOpt = res$lambda[indLambdaOpt], selectedGroups = REJECT[[indLambdaOpt[1]]], lambda = res$lambda, nGroup = res$nGroup, nSelectedGroup = nbReject, 
              group = group, var = var, test = test, alpha = alpha, reject = REJECT, control = control, time = time, hierTest = TEST)
  
  class(out) = "HMT"
  
  return(out)
}