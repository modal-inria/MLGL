#'
#' Obtain a sparse matrix of the coefficients of the path
#' 
#' @param x \code{\link{MLGL}} object
#' @param row "lambda" or "covariates". If row="covariates", each row of the output matrix represents a covariates else ff row="lambda", it represents a value of lambda.
#'
#' @return a sparse matrix containing the estimated coefficients for different values of lambda
#'
#' @details This functions can be used with a \code{\link{MLGL}} object to obtain a matrix with all estimated coefficients for the p original variables.
#' In case of overlapping groups, coefficients from repeated variables are summed. 
#'
#' @examples 
#' # Simulate gaussian data with block-diagonal variance matrix containing 12 blocks of size 5
#' X <- simuBlockGaussian(50, 12, 5, 0.7)
#' # Generate a response variable
#' y <- X[,c(2,7,12)]%*%c(2,2,-2) + rnorm(50,0,0.5)
#' # Apply MLGL method
#' res <- MLGL(X,y)
#' # Convert output in sparse matrix format
#' beta <- listToMatrix(res)
#'
#' @seealso \link{MLGL}, \link{overlapgglasso}
#'
#' @export
listToMatrix <- function(x, row = c("covariates","lambda"))
{
  row = match.arg(row)
  if(row == "covariates")
  {
    bet <- Matrix(0, ncol = length(x$lambda), nrow = x$dim[2])
    for(i in 1:length(x$lambda))
    {
      dup <- duplicated(x$var[[i]])
      bet[x$var[[i]][!dup], i] = bet[x$var[[i]][!dup],i] + x$beta[[i]][!dup]
      bet[x$var[[i]][dup] , i] = bet[x$var[[i]][dup] ,i] + x$beta[[i]][dup]
    }
    
    return(bet)
  }
  else
  {
    bet <- Matrix(0, nrow = length(x$lambda), ncol = x$dim[2])
    for(i in 1:length(x$lambda))
    {
      dup <- duplicated(x$var[[i]])
      bet[i,x$var[[i]][!dup]] = bet[i,x$var[[i]][!dup]] + x$beta[[i]][!dup]
      bet[i,x$var[[i]][dup]]  = bet[i,x$var[[i]][dup]]  + x$beta[[i]][dup]
    }
    
    return(bet)
  }

}
