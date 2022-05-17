#'
#' Predict fitted values from a \code{\link{MLGL}} object
#'
#' @author original code from \pkg{gglasso} package Author: Yi Yang <yiyang@@umn.edu>, Hui Zou <hzou@@stat.umn.edu>
#'
#' @param object \code{\link{MLGL}} object
#' @param newx matrix with new individuals for prediction. If type="coefficients", the parameter has to be NULL
#' @param s values of lambda. If NULL, use values from object
#' @param type if "fit", return the fitted values for each values of s, if "coefficients", return the estimated 
#' coefficients for each s
#' @param ... Not used. Other arguments to predict.
#'
#' @return A matrix with fitted values or estimated coefficients for given values of s.
#'
#' @method predict MLGL
#'
#' @author function inspired from predict function from gglasso package by Yi Yang and Hui Zou.
#'
#' @examples 
#' X <- simuBlockGaussian(n = 50, nBlock = 12, sizeBlock = 5, rho = 0.7)
#' y <- drop(X[, c(2, 7, 12)] %*% c(2, 2, -1)) + rnorm(50, 0, 0.5)
#' 
#' m1 <- MLGL(X, y, loss = "ls")
#' predict(m1, newx = X)
#' predict(m1, s=3, newx = X)
#' predict(m1, s=1:3, newx = X)
#'
#' @seealso \link{MLGL}
#'
#' @export
predict.MLGL <- function(object, newx = NULL, s = NULL, type = c("fit", "coefficients"), ...) {
  # check values of type
  type <- match.arg(type)

  if ((type == "fit") & is.null(newx)) {
    stop("newx is missing with type='fit'.")
  }
  
  if (!is.null(s)) {
    if (!is.numeric(s)) {
      stop("s must be a vector of positive real.")
    }
    if (any(s < 0)) {
      stop("s must be a vector of positive real.")
    }
  }

  b0 <- t(as.matrix(object$b0))
  rownames(b0) <- "(Intercept)"
  nbeta <- rbind2(b0, listToMatrix(object))
  if (!is.null(s)) {
    vnames <- dimnames(nbeta)[[1]]
    dimnames(nbeta) <- list(NULL, NULL)
    lambda <- object$lambda
    lamlist <- lambda.interp(lambda, s)
    if (length(s) == 1) {
      nbeta <- nbeta[, lamlist$left, drop = FALSE] * lamlist$frac + nbeta[, lamlist$right, drop = FALSE] * (1 - lamlist$frac)
    }
    else {
      nbeta <- nbeta[, lamlist$left, drop = FALSE] %*% diag(lamlist$frac) + nbeta[, lamlist$right, drop = FALSE] %*% diag(1 - lamlist$frac)
    }
    dimnames(nbeta) <- list(vnames, paste(seq(along = s)))
  }

  if (!is.null(newx) & is.null(dim(newx))) {
    newx <- matrix(newx, nrow = 1)
  }

  return(switch(type, fit = as.matrix(as.matrix(cbind2(1, newx)) %*% nbeta), coefficients = nbeta))
}

#'
#' Predict fitted values from a \code{\link{cv.MLGL}} object
#'
#' @author Quentin Grimonprez
#'
#' @param object \code{\link{cv.MLGL}} object
#' @param newx matrix with new individuals for prediction. If type="coefficients", the parameter has to be NULL
#' @param s Either "lambda.1se" or "lambda.min"
#' @param type if "fit", return the fitted values for each values of s, if "coefficients", return the estimated 
#' coefficients for each s
#' @param ... Not used. Other arguments to predict.
#'
#' @return A matrix with fitted values or estimated coefficients for given values of s.
#'
#' @method predict cv.MLGL
#'
#' @seealso \link{cv.MLGL}
#'
#' @export
predict.cv.MLGL <- function(object, newx = NULL, s = c("lambda.1se", "lambda.min"), type = c("fit", "coefficients"), ...) {
  # check s
  s <- match.arg(s)

  # check type
  type <- match.arg(type)

  # get the desired value of lambda
  lambda <- switch(s,
    lambda.1se = object$lambda.1se,
    lambda.min = object$lambda.min
  )

  return(predict(object, newx, s = lambda, type = "fit"), ...)
}


#'
#' Get coefficients from a \code{\link{MLGL}} object
#'
#' @author Quentin Grimonprez
#'
#' @param object \code{\link{MLGL}} object
#' @param s values of lambda. If NULL, use values from object
#' @param ... Not used. Other arguments to predict.
#'
#' @return A matrix with estimated coefficients for given values of s.
#'
#' @method coef MLGL
#'
#' @seealso \link{MLGL}, \link{predict.MLGL}
#'
#' @export
coef.MLGL <- function(object, s = NULL, ...) {
  return(predict(object, s = s, type = "coefficients", ...))
}


#'
#' Get coefficients from a \code{\link{cv.MLGL}} object
#'
#' @author Quentin Grimonprez
#'
#' @param object \code{\link{cv.MLGL}} object
#' @param s Either "lambda.1se" or "lambda.min"
#' @param ... Not used. Other arguments to predict.
#'
#' @return A matrix with estimated coefficients for given values of s.
#'
#' @method coef cv.MLGL
#'
#' @seealso \link{cv.MLGL}, \link{predict.cv.MLGL}
#'
#' @export
coef.cv.MLGL <- function(object, s = c("lambda.1se", "lambda.min"), ...) {
  # check s
  s <- match.arg(s)

  # get the desired value of lambda
  lambda <- switch(s, lambda.1se = object$lambda.1se, lambda.min = object$lambda.min)

  return(predict(object, s = lambda, type = "coefficients"))
}



### code from gglasso package
### Author: Yi Yang <yiyang@umn.edu>, Hui Zou <hzou@stat.umn.edu>
lambda.interp <- function(lambda, s) {
  ### lambda is the index sequence that is produced by the model
  ### s is the new vector at which evaluations are required.
  ### the value is a vector of left and right indices, and a
  #   vector of fractions.
  ### the new values are interpolated between the two using the
  #   fraction
  ### Note: lambda decreases. you take:
  ### sfrac*left+(1-sfrac*right)
  if (length(lambda) == 1) {
    left <- rep(1, length(s))
    right <- left
    sfrac <- rep(1, length(s))
  } else {
    s[s > max(lambda)] <- max(lambda)
    s[s < min(lambda)] <- min(lambda)
    k <- length(lambda)
    sfrac <- (lambda[1] - s) / (lambda[1] - lambda[k])
    lambda <- (lambda[1] - lambda) / (lambda[1] - lambda[k])
    coord <- approx(lambda, seq(lambda), sfrac)$y
    left <- floor(coord)
    right <- ceiling(coord)
    sfrac <- (sfrac - lambda[right]) / (lambda[left] - lambda[right])
    sfrac[left == right] <- 1
  }
  list(left = left, right = right, frac = sfrac)
}
