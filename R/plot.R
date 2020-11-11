
#'
#' Plot the path obtained from \code{\link{MLGL}} function
#'
#' @param x \code{\link{MLGL}} object
#' @param log.lambda If TRUE, use log(lambda) instead of lambda in abscissa
#' @param lambda.lines if TRUE, add vertical lines at lambda values
#' @param ... Other parameters for plot function
#'
#' @examples
#' # Simulate gaussian data with block-diagonal variance matrix containing 12 blocks of size 5
#' set.seed(42)
#' X <- simuBlockGaussian(50, 12, 5, 0.7)
#' # Generate a response variable
#' y <- X[, c(2, 7, 12)] %*% c(2, 2, -2) + rnorm(50, 0, 0.5)
#' # Apply MLGL method
#' res <- MLGL(X, y)
#' # Plot the solution path
#' plot(res)
#' @method plot MLGL
#'
#' @seealso \link{MLGL}
#'
#' @export
plot.MLGL <- function(x, log.lambda = FALSE, lambda.lines = FALSE, ...) {
  # check log
  if (length(log.lambda) != 1) {
    stop("log must be a boolean.")
  }
  if (!is.logical(log.lambda)) {
    stop("log must be a boolean.")
  }
  # check log
  if (length(lambda.lines) != 1) {
    stop("lambda.lines must be a boolean.")
  }
  if (!is.logical(lambda.lines)) {
    stop("lambda.lines must be a boolean.")
  }

  bet <- listToMatrix(x, "lambda")

  # abscissa : log or not ?
  absc <- x$lambda
  if (log.lambda) {
    absc <- log(absc)
  }

  # plot
  matplot(absc, bet,
    type = "l", lty = 1, xlab = ifelse(log.lambda, expression(paste("log(", lambda, ")")), expression(lambda)),
    ylab = "Coefficients", main = "Solution path", ...
  )

  # add vertical lines for lambda values
  if (lambda.lines) {
    abline(v = absc, col = "blue", lty = "dotted", lwd = 0.5)
  }
}


#'
#' Plot the cross-validation obtained from \code{\link{cv.MLGL}} function
#'
#' @param x \code{\link{cv.MLGL}} object
#' @param log.lambda If TRUE, use log(lambda) instead of lambda in abscissa
#' @param ... Other parameters for plot function
#'
#' @examples
#' set.seed(42)
#' # Simulate gaussian data with block-diagonal variance matrix containing 12 blocks of size 5
#' X <- simuBlockGaussian(50, 12, 5, 0.7)
#' # Generate a response variable
#' y <- X[, c(2, 7, 12)] %*% c(2, 2, -2) + rnorm(50, 0, 0.5)
#' # Apply cv.MLGL method
#' res <- cv.MLGL(X, y)
#' # Plot the cv error curve
#' plot(res)
#' @method plot cv.MLGL
#'
#' @seealso \link{cv.MLGL}
#'
#' @export
plot.cv.MLGL <- function(x, log.lambda = FALSE, ...) {
  # check log
  if (length(log.lambda) != 1) {
    stop("log must be a boolean.")
  }
  if (!is.logical(log.lambda)) {
    stop("log must be a boolean.")
  }

  # abscissa : log or not ?
  absc <- x$lambda
  lam <- c(x$lambda.min, x$lambda.1se)
  if (log.lambda) {
    absc <- log(absc)
    lam <- log(lam)
  }

  # plot
  matplot(absc, cbind(x$cvm, x$cvupper, x$cvlower),
    type = "l", lty = c(1, 2, 2), col = c(1, 2, 2),
    xlab = ifelse(log.lambda, expression(paste("log(", lambda, ")")), expression(lambda)), ylab = "Error", ...
  )
  abline(v = lam, col = "blue", lty = "dashed")
}


#'
#' Plot the stability path obtained from \code{\link{stability.MLGL}} function
#'
#' @param x \code{\link{stability.MLGL}} object
#' @param log.lambda If TRUE, use log(lambda) instead of lambda in abscissa
#' @param threshold Threshold for selection frequency
#' @param ... Other parameters for plot function
#'
#' @return A list containing :
#' \describe{
#' \item{var}{Index of selected variables for the given threshold.}
#' \item{group}{Index of the associated group.}
#' \item{threshold}{Value of threshold}
#' }
#'
#' @examples
#' \donttest{
#' set.seed(42)
#' # Simulate gaussian data with block-diagonal variance matrix containing 12 blocks of size 5
#' X <- simuBlockGaussian(50, 12, 5, 0.7)
#'
#' # Generate a response variable
#' y <- X[, c(2, 7, 12)] %*% c(2, 2, -2) + rnorm(50, 0, 0.5)
#'
#' # Apply stability.MLGL method
#' res <- stability.MLGL(X, y)
#' selected <- plot(res)
#' print(selected)
#' }
#'
#' @method plot stability.MLGL
#'
#' @seealso \link{stability.MLGL}
#'
#' @export
plot.stability.MLGL <- function(x, log.lambda = FALSE, threshold = 0.75, ...) {
  # check log
  if (length(log.lambda) != 1) {
    stop("log must be a boolean.")
  }
  if (!is.logical(log.lambda)) {
    stop("log must be a boolean.")
  }
  # threshold
  if (!is.numeric(threshold) | (length(threshold) != 1)) {
    stop("threshold must be a positive real lesser than 1.")
  }
  if ((threshold < 0) | (threshold > 1)) {
    stop("threshold must be a positive real lesser than 1.")
  }

  # abscissa : log or not ?
  absc <- x$lambda
  if (log.lambda) {
    absc <- log(absc)
  }

  # determine color according to threshold
  col <- apply(x$stability, 2, FUN = function(x) {
    ifelse(any(x > threshold), 2, 1)
  })

  # plot
  matplot(absc, x$stability,
    type = "l", lty = 1, col = col,
    xlab = ifelse(log.lambda, expression(paste("log(", lambda, ")")), expression(lambda)), ylab = "Probability selection", ...
  )
  abline(h = threshold, col = "blue", lty = "dashed")

  # determine selected groups and variables
  selectedGroup <- which(col == 2)
  indsel <- x$group %in% selectedGroup

  return(list(var = x$var[indsel], group = x$group[indsel], threshold = threshold))
}



#'
#' Plot the path obtained from \code{\link{fullProcess}} function
#'
#' @param x \code{\link{fullProcess}} object
#' @param log.lambda If TRUE, use log(lambda) instead of lambda in abscissa
#' @param lambda.lines If TRUE, add vertical lines at lambda values
#' @param lambda.opt If there is several optimal lambdas, which one to print "min", "max" or "both"
#' @param ... Other parameters for plot function
#'
#' @examples
#' set.seed(42)
#' # Simulate gaussian data with block-diagonal variance matrix containing 12 blocks of size 5
#' X <- simuBlockGaussian(50, 12, 5, 0.7)
#' # Generate a response variable
#' y <- X[, c(2, 7, 12)] %*% c(2, 2, -2) + rnorm(50, 0, 0.5)
#' # Apply MLGL method
#' res <- fullProcess(X, y)
#' # Plot the solution path
#' plot(res)
#' @method plot fullProcess
#'
#' @seealso \link{fullProcess}
#'
#' @export
plot.fullProcess <- function(x, log.lambda = FALSE, lambda.lines = FALSE, lambda.opt = c("min", "max", "both"), ...) {
  lambda.opt <- match.arg(lambda.opt)

  lambdaOpt <- switch(lambda.opt,
    "min" = x$lambdaOpt[1],
    "max" = x$lambdaOpt[length(x$lambdaOpt)],
    "both" = c(1, x$lambdaOpt[length(x$lambdaOpt)])
  )
  par(mfrow = c(2, 1))
  plot(x$res, ...)
  abline(v = ifelse(log.lambda, log(lambdaOpt), lambdaOpt), col = "red", lty = "dotted")
  text(ifelse(log.lambda, log(lambdaOpt), lambdaOpt), 0, labels = expression(lambda[opt]), col = "red", pos = 1)

  abscissa <- x$res$lambda
  if (log.lambda) {
    abscissa <- log(abscissa)
  }

  matplot(abscissa, cbind(x$res$nGroup, sapply(x$reject, length)),
    type = "l", col = c(1, 4), xlab = ifelse(log.lambda, expression(paste("log(", lambda, ")")), expression(lambda)),
    ylab = "Number of groups", main = "Number of groups in MLGL path", ...
  )
  abline(v = ifelse(log.lambda, log(lambdaOpt), lambdaOpt), col = "red", lty = "dotted")
  text(ifelse(log.lambda, log(lambdaOpt), lambdaOpt), 0, labels = expression(lambda[opt]), col = "red", pos = 2)
  legend("topright", c("before testing", "after testing"), col = c(1, 4), lty = 1:2, cex = 0.6)
}


#'
#' Plot the path obtained from \code{\link{HMT}} function
#'
#' @param x \code{\link{fullProcess}} object
#' @param log.lambda If TRUE, use log(lambda) instead of lambda in abscissa
#' @param lambda.lines If TRUE, add vertical lines at lambda values
#' @param lambda.opt If there is several optimal lambdas, which one to print "min", "max" or "both"
#' @param ... Other parameters for plot function
#'
#' @examples
#' set.seed(42)
#' # Simulate gaussian data with block-diagonal variance matrix containing 12 blocks of size 5
#' X <- simuBlockGaussian(50, 12, 5, 0.7)
#' # Generate a response variable
#' y <- X[, c(2, 7, 12)] %*% c(2, 2, -2) + rnorm(50, 0, 0.5)
#' # Apply MLGL method
#' res <- MLGL(X, y)
#'
#' out <- HMT(res, X, y)
#' plot(out)
#' @method plot HMT
#'
#' @seealso \link{HMT}
#'
#' @export
plot.HMT <- function(x, log.lambda = FALSE, lambda.lines = FALSE, lambda.opt = c("min", "max", "both"), ...) {
  lambda.opt <- match.arg(lambda.opt)

  lambdaOpt <- switch(lambda.opt,
    "min" = x$lambdaOpt[1],
    "max" = x$lambdaOpt[length(x$lambdaOpt)],
    "both" = c(1, x$lambdaOpt[length(x$lambdaOpt)])
  )

  abscissa <- x$lambda
  if (log.lambda) {
    abscissa <- log(abscissa)
  }

  matplot(abscissa, cbind(x$nGroup, x$nSelectedGroup),
    type = "l", col = c(1, 4), xlab = ifelse(log.lambda, expression(paste("log(", lambda, ")")), expression(lambda)),
    ylab = "Number of groups", main = "Number of groups in MLGL path", ...
  )
  abline(v = ifelse(log.lambda, log(lambdaOpt), lambdaOpt), col = "red", lty = "dotted")
  text(ifelse(log.lambda, log(lambdaOpt), lambdaOpt), 0, labels = expression(lambda[opt]), col = "red", pos = 2)
  legend("topright", c("before testing", "after testing"), col = c(1, 4), lty = 1:2, cex = 0.6)
}
