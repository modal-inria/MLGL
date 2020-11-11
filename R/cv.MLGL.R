#' V-fold cross validation for \code{\link{MLGL}} function
#'
#' @title  Multi-Layer Group-Lasso with cross V-fold validation
#' @author Quentin Grimonprez
#' @param X matrix of size n*p
#' @param y vector of size n. If loss = "logit", elements of y must be in {-1,1}
#' @param nfolds number of folds
#' @param lambda lambda values for group lasso. If not provided, the function generates its own values of lambda
#' @param hc output of \code{\link{hclust}} function. If not provided, \code{\link{hclust}} is run with ward.D2 method
#' @param weightLevel a vector of size p for each level of the hierarchy. A zero indicates that the level will be ignored. If not provided, use 1/(height between 2 successive levels)
#' @param weightSizeGroup a vector
#' @param loss a character string specifying the loss function to use, valid options are: "ls" least squares loss (regression) and "logit" logistic loss (classification)
#' @param sizeMaxGroup maximum size of selected groups. If NULL, no restriction
#' @param intercept should an intercept be included in the model ?
#' @param verbose print some informations
#' @param ... Others parameters for \code{\link{cv.gglasso}} function
#'
#'
#' @return a cv.MLGL object containing :
#' \describe{
#' \item{lambda}{values of \code{lambda}.}
#' \item{cvm}{the mean cross-validated error.}
#' \item{cvsd}{estimate of standard error of \code{cvm}}
#' \item{cvupper}{upper curve = \code{cvm+cvsd}}
#' \item{cvlower}{lower curve = \code{cvm-cvsd}}
#' \item{lambda.min}{The optimal value of \code{lambda} that gives minimum cross validation error \code{cvm}.}
#' \item{lambda.1se}{The largest value of \code{lambda} such that error is within 1 standard error of the minimum.}
#' \item{time}{computation time}
#' }
#'
#'
#' @details
#' Hierarhical clustering is performed with all the variables. Then, the partitions from the different
#'  levels of the hierarchy are used in the differents run of MLGL for cross validation.
#'
#' @examples
#' set.seed(42)
#' # Simulate gaussian data with block-diagonal variance matrix containing 12 blocks of size 5
#' X <- simuBlockGaussian(50, 12, 5, 0.7)
#' # Generate a response variable
#' y <- X[, c(2, 7, 12)] %*% c(2, 2, -2) + rnorm(50, 0, 0.5)
#' # Apply cv.MLGL method
#' res <- cv.MLGL(X, y)
#' @seealso \link{MLGL}, \link{stability.MLGL}, \link{predict.cv.gglasso}, \link{coef.cv.MLGL}, \link{plot.cv.MLGL}
#'
#' @export
cv.MLGL <- function(X, y, nfolds = 5, lambda = NULL, hc = NULL, weightLevel = NULL, weightSizeGroup = NULL, loss = c("ls", "logit"), intercept = TRUE, sizeMaxGroup = NULL, verbose = FALSE, ...) {

  # check parameters
  loss <- match.arg(loss)
  .checkParameters(X, y, hc, lambda, weightLevel, weightSizeGroup, intercept, verbose, loss, sizeMaxGroup)

  # nfolds
  if (!is.numeric(nfolds) | (length(nfolds) != 1)) {
    stop("nfolds must be an integer greater than 3.")
  }
  if (!.is.wholenumber(nfolds) | nfolds <= 2 | nfolds > nrow(X)) { # restriction >=3 comes from cv.gglasso function
    stop("nfolds must be an integer greater than 3.")
  }

  ################################ same as MLGL function
  # define some usefull variables
  n <- nrow(X)
  p <- ncol(X)
  tcah <- rep(NA, 3)

  # if no hc output provided, we make one
  if (is.null(hc)) {
    if (verbose) {
      cat("Computing hierarchical clustering...")
    }

    t1 <- proc.time()
    d <- dist(t(X))
    hc <- hclust(d, method = "ward.D2")
    t2 <- proc.time()
    tcah <- t2 - t1
    if (verbose) {
      cat("DONE in ", tcah[3], "s\n")
    }
  }

  # compute weight, active variables and groups
  if (verbose) {
    cat("Preliminary step...")
  }
  t1 <- proc.time()
  prelim <- preliminaryStep(hc, weightLevel, weightSizeGroup, sizeMaxGroup)

  # duplicate data
  Xb <- X[, prelim$var]
  t2 <- proc.time()
  if (verbose) {
    cat("DONE in ", (t2 - t1)[3], "s\n")
  }

  ################################ END same as MLGL function



  ######## group lasso
  if (verbose) {
    cat("Computing group-lasso...")
  }
  t1 <- proc.time()
  res <- cv.gglasso(Xb, y, prelim$group, pf = prelim$weight, nfolds = nfolds, lambda = lambda, intercept = intercept, loss = loss, ...)
  t2 <- proc.time()
  tgglasso <- t2 - t1
  if (verbose) {
    cat("DONE in ", tgglasso[3], "s\n")
  }

  # delete some gglasso output
  res$name <- NULL
  res$gglasso.fit <- NULL

  # rename output
  res$cvlower <- res$cvlo
  res$cvlo <- NULL

  res$time <- c(tcah[3], tgglasso[3])
  names(res$time) <- c("hclust", "glasso")
  class(res) <- "cv.MLGL"


  return(res)
}
