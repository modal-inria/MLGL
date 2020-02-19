
#' simulate n samples from a gaussian multivariate law with 0 vector mean and block diagonal variance matrix with diagonal 1 and block of rho.
#'
#' @title Simulate multivariate Gaussian samples with block diagonal varaince matrix
#'
#' @author Quentin Grimonprez
#' @param n number of samples to simulate
#' @param nBlock number of blocks
#' @param sizeBlock size of blocks
#' @param rho correlation within each block
#'
#' @return a matrix of size n *(nBlock*sizeBlock) containing the samples
#'
#' @examples 
#' X = simuBlockGaussian(50,12,5,0.7)
#' @export
simuBlockGaussian <- function(n, nBlock, sizeBlock, rho)
{
  .checkSimu(n,nBlock,sizeBlock,rho)
  
  Sigma <- matrix(rho, nrow = sizeBlock, ncol = sizeBlock)
  diag(Sigma)=1
  X <- matrix(nrow = n, ncol = nBlock*sizeBlock)
  
  for(i in 1:nBlock)
    X[,((i-1)*sizeBlock+1):(i*sizeBlock)] = MASS::mvrnorm(n, rep(0,sizeBlock), Sigma)
  
  return(X)
}

#check parameters from simuBlockGaussian function
.checkSimu <- function(n, nBlock, sizeBlock, rho)
{
  #number of samples
  if( !is.numeric(n) | (length(n)!=1) )
    stop("n must be a strictly positive integer.")
  if(!.is.wholenumber(n) | n<=0)
    stop("n must be a strictly positive integer.")
  
  #nBlock
  if( !is.numeric(nBlock) | (length(nBlock)!=1) )
    stop("nBlock must be a strictly positive integer.")
  if(!.is.wholenumber(nBlock) | nBlock<=0)
    stop("nBlock must be a strictly positive integer.")
  
  #sizeBlock
  if( !is.numeric(sizeBlock) | (length(sizeBlock)!=1) )
    stop("sizeBlock must be a strictly positive integer.")
  if(!.is.wholenumber(sizeBlock) | sizeBlock<=0)
    stop("sizeBlock must be a strictly positive integer.")
  
  #rho
  if( !is.numeric(rho) | (length(rho)!=1) )
    stop("rho must be a real between -1 and 1.")
  if(abs(rho)>1)
    stop("rho must be a real between -1 and 1.")
}

#check if a number is an integer
.is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  
{
  abs(x - round(x)) < tol
}