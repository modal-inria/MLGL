#' @title Print Values
#'
#' Print a \code{\link{MLGL}} object
#'
#'
#' @param x \code{\link{MLGL}} object
#' @param ... Not used.
#'
#'
#' @method print MLGL
#'
#' @examples
#' set.seed(42)
#' # Simulate gaussian data with block-diagonal variance matrix containing 12 blocks of size 5
#' X <- simuBlockGaussian(50, 12, 5, 0.7)
#' # Generate a response variable
#' y <- X[, c(2, 7, 12)] %*% c(2, 2, -2) + rnorm(50, 0, 0.5)
#' # Apply MLGL method
#' res <- MLGL(X, y)
#' print(res)
#' @seealso \link{MLGL} \link{summary.MLGL}
#'
#' @export
print.MLGL <- function(x, ...) {
  cat("$lambda\n")
  print(x$lambda)
  cat("$nVar\n")
  print(x$nVar)
  cat("$nGroup\n")
  print(x$nGroup)
}

#' @title Object Summaries
#'
#' Summary of a \code{\link{MLGL}} object
#'
#'
#' @param object \code{\link{MLGL}} object
#' @param ... Not used.
#'
#'
#' @method summary MLGL
#'
#' @examples
#' set.seed(42)
#' # Simulate gaussian data with block-diagonal variance matrix containing 12 blocks of size 5
#' X <- simuBlockGaussian(50, 12, 5, 0.7)
#' # Generate a response variable
#' y <- X[, c(2, 7, 12)] %*% c(2, 2, -2) + rnorm(50, 0, 0.5)
#' # Apply MLGL method
#' res <- MLGL(X, y)
#' summary(res)
#' @seealso \link{MLGL} \link{print.MLGL}
#'
#' @export
summary.MLGL <- function(object, ...) {
  cat("#### MLGL\n")
  cat("## Data \n")
  cat("Number of individuals:", object$dim[1], "\n")
  cat("Number of variables:", object$dim[2], "\n")
  cat("\n")
  cat("## Hierarchical clustering \n")
  cat("HC proveded by user:", "hc" %in% names(object$call), "\n")
  cat("Time:", object$time[1], "s\n")
  cat("\n")
  cat("## Group-lasso\n")
  cat("Loss:", object$loss, "\n")
  cat("Intercept:", object$intercept, "\n")
  cat("Number of lambda:", length(object$lambda), "\n")
  cat("Number of selected variables:", head(object$nVar), "...\n")
  cat("Number of selected groups:", head(object$nGroup), "...\n")
  cat("Time:", object$time[2], "s\n")
  cat("\n")
  cat("Total elapsed time:", sum(object$time, na.rm = TRUE), "s\n")
}

#' @title Print Values
#'
#' Print a \code{\link{fullProcess}} object
#'
#'
#' @param x \code{\link{fullProcess}} object
#' @param ... Not used.
#'
#' @examples
#' set.seed(42)
#' # Simulate gaussian data with block-diagonal variance matrix containing 12 blocks of size 5
#' X <- simuBlockGaussian(50, 12, 5, 0.7)
#' # Generate a response variable
#' y <- X[, c(2, 7, 12)] %*% c(2, 2, -2) + rnorm(50, 0, 0.5)
#' # Apply MLGL method
#' res <- fullProcess(X, y)
#' print(res)
#' @method print fullProcess
#'
#' @seealso \link{fullProcess} \link{summary.fullProcess}
#'
#' @export
print.fullProcess <- function(x, ...) {
  cat("Group-lasso\n")
  cat("$res$lambda\n")
  print(x$res$lambda)
  cat("$res$nVar\n")
  print(x$res$nVar)
  cat("$res$nGroup\n")
  print(x$res$nGroup)
  cat("Test output\n")
  cat("$lambdaOpt\n")
  print(x$lambdaOpt)
  cat("$selectedGroups\n")
  cat(x$selectedGroups)
}

#' @title Object Summaries
#'
#' Summary of a \code{\link{fullProcess}} object
#'
#'
#' @param object \code{\link{fullProcess}} object
#' @param ... Not used.
#'
#' @examples
#' set.seed(42)
#' # Simulate gaussian data with block-diagonal variance matrix containing 12 blocks of size 5
#' X <- simuBlockGaussian(50, 12, 5, 0.7)
#' # Generate a response variable
#' y <- X[, c(2, 7, 12)] %*% c(2, 2, -2) + rnorm(50, 0, 0.5)
#' # Apply MLGL method
#' res <- fullProcess(X, y)
#' summary(res)
#' @method summary fullProcess
#'
#' @seealso \link{fullProcess} \link{print.fullProcess}
#'
#' @export
summary.fullProcess <- function(object, ...) {
  summary(object$res)
  cat("#### Multiple Hierarchical testing\n")
  cat("## Data \n")
  cat("alpha:", object$alpha, "\n")
  cat("control:", object$control, "\n")
  cat("optimal lambda:\n")
  print(object$lambdaOpt)
  cat("Selected groups:", object$selectedGroups, "\n")
  cat("Selected variables:\n")
  print(object$var)
  cat("Time:", object$time[3], "s\n")
  cat("\n")
  cat("Total elapsed time:", sum(object$time, na.rm = TRUE), "s\n")
}



#' @title Print Values
#'
#' Print a \code{\link{HMT}} object
#'
#'
#' @param x \code{\link{HMT}} object
#' @param ... Not used.
#'
#'
#' @examples
#' set.seed(42)
#' # Simulate gaussian data with block-diagonal variance matrix containing 12 blocks of size 5
#' X <- simuBlockGaussian(50, 12, 5, 0.7)
#' # Generate a response variable
#' y <- X[, c(2, 7, 12)] %*% c(2, 2, -2) + rnorm(50, 0, 0.5)
#' # Apply MLGL method
#' res <- MLGL(X, y)
#' out <- HMT(res, X, y)
#' print(out)
#' @method print HMT
#'
#' @seealso \link{HMT} \link{summary.HMT}
#'
#' @export
print.HMT <- function(x, ...) {
  cat("$lambda\n")
  print(x$lambda)
  cat("$nGroup\n")
  print(x$nGroup)
  cat("Test output\n")
  cat("$nSelectedGroup\n")
  print(x$nSelectedGroup)
  cat("$lambdaOpt\n")
  cat(x$lambdaOpt)
  cat("$selectedGroups\n")
  cat(x$selectedGroups)
}

#' @title Object Summaries
#'
#' Summary of a \code{\link{HMT}} object
#'
#'
#' @param object \code{\link{HMT}} object
#' @param ... Not used.
#'
#' @examples
#' set.seed(42)
#' # Simulate gaussian data with block-diagonal variance matrix containing 12 blocks of size 5
#' X <- simuBlockGaussian(50, 12, 5, 0.7)
#' # Generate a response variable
#' y <- X[, c(2, 7, 12)] %*% c(2, 2, -2) + rnorm(50, 0, 0.5)
#' # Apply MLGL method
#' res <- MLGL(X, y)
#' out <- HMT(res, X, y)
#' summary(out)
#' @method summary HMT
#'
#' @seealso \link{HMT} \link{print.HMT}
#'
#' @export
summary.HMT <- function(object, ...) {
  cat("#### Multiple Hierarchical testing\n")
  cat("## Data \n")
  cat("alpha:", object$alpha, "\n")
  cat("control:", object$control, "\n")
  cat("optimal lambda:", object$lambdaOpt, "\n")
  cat("Selected groups:", object$selectedGroups, "\n")
  cat("Selected variables:", object$var, "\n")
  cat("Time:", object$time[3], "s\n")
  cat("\n")
  cat("Total elapsed time:", sum(object$time, na.rm = TRUE), "s\n")
}
