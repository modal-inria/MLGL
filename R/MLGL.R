#' Run hierarchical clustering following by a group-lasso on all the different partitions.
#'
#' @title Multi-Layer Group-Lasso
#'
#' @author Quentin Grimonprez
#' @param X matrix of size n*p
#' @param y vector of size n. If loss = "logit", elements of y must be in {-1,1}
#' @param hc output of \code{\link{hclust}} function. If not provided, \code{\link{hclust}} is run with \code{ward.D2} method. User can also provide the desired method: "single", "complete", "average", "mcquitty", "ward.D", "ward.D2", "centroid", "median".
#' @param lambda lambda values for group lasso. If not provided, the function generates its own values of lambda
#' @param weightLevel a vector of size p for each level of the hierarchy. A zero indicates that the level will be ignored. If not provided, use 1/(height between 2 successive levels). Only if \code{hc} is provided
#' @param weightSizeGroup a vector of size 2*p-1 containing the weight for each group. Default is the square root of the size of each group. Only if \code{hc} is provided
#' @param intercept should an intercept be included in the model ?
#' @param loss a character string specifying the loss function to use, valid options are: "ls" least squares loss (regression) and "logit" logistic loss (classification)
#' @param sizeMaxGroup maximum size of selected groups. If NULL, no restriction
#' @param verbose print some information
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
#'   weight : weight associated with the different groups.
#'   level : for each group, the corresponding level of the hierarchy where it appears and disappears. 3 indicates the level with a partition of 3 groups.}
#'   \item{time}{computation time}
#'   \item{dim}{dimension of \code{X}}
#'   \item{hc}{Output of hierarchical clustering}
#'   \item{call}{Code executed by user}
#' }
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
#' @seealso \link{cv.MLGL}, \link{stability.MLGL}, \link{listToMatrix}, \link{predict.MLGL}, \link{coef.MLGL}, \link{plot.cv.MLGL}
#'
#' @export
MLGL <- function(X, y, hc = NULL, lambda = NULL, weightLevel = NULL, weightSizeGroup = NULL, intercept = TRUE, loss = c("ls", "logit"), sizeMaxGroup = NULL, verbose = FALSE, ...) {
  # check parameters
  loss <- match.arg(loss)
  .checkParameters(X, y, hc, lambda, weightLevel, weightSizeGroup, intercept, verbose, loss, sizeMaxGroup)

  # define some usefull variables
  n <- nrow(X)
  p <- ncol(X)
  tcah <- NA

  ######## hierarchical clustering
  # if hc output not provided, we perform one
  if (is.null(hc) | is.character(hc)) {
    if (verbose) {
      cat("Computing hierarchical clustering...")
    }

    t1 <- proc.time()
    d <- dist(t(X))
    hc <- fastcluster::hclust(d, method = ifelse(is.character(hc), hc, "ward.D2"))
    t2 <- proc.time()

    hc$time <- as.numeric((t2 - t1)[3])
    tcah <- hc$time

    if (verbose) {
      cat("DONE in ", tcah, "s\n")
    }
  }

  ######## compute weight, active variables and groups
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

  ######## group lasso
  if (verbose) {
    cat("Computing group-lasso...")
  }
  t1 <- proc.time()
  res <- gglasso(Xb, y, prelim$group, pf = prelim$weight, lambda = lambda, intercept = intercept, loss = loss, ...)
  t2 <- proc.time()
  tgglasso <- as.numeric((t2 - t1)[3])
  if (verbose) {
    cat("DONE in ", tgglasso, "s\n")
  }

  ########  create output object
  res2 <- list()
  res2$lambda <- res$lambda
  non0 <- apply(res$beta, 2, FUN = function(x) {
    which(x != 0)
  })
  res2$var <- lapply(non0, FUN = function(x) {
    prelim$var[x]
  })
  res2$nVar <- sapply(res2$var, FUN = function(x) {
    length(unique(x))
  })
  res2$group <- lapply(non0, FUN = function(x) {
    prelim$group[x]
  })
  res2$nGroup <- sapply(res2$group, FUN = function(x) {
    length(unique(x))
  })
  res2$beta <- lapply(1:length(res$lambda), FUN = function(x) {
    res$beta[non0[[x]], x]
  })
  res2$b0 <- res$b0
  res2$structure <- prelim
  res2$dim <- dim(X)
  res2$hc <- hc
  res2$time <- c(tcah, tgglasso)
  names(res2$time) <- c("hclust", "glasso")
  res2$call <- match.call()
  res2$intercept <- intercept
  res2$loss <- loss
  class(res2) <- "MLGL"


  return(res2)
}



#' @param formula an object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted.
#' @param data an optional data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model. If not found in data, the variables are taken from environment(formula)
#'
#'
#' @rdname MLGL
#'
#' @export
MLGL.formula <- function(formula, data, hc = NULL, lambda = NULL, weightLevel = NULL, weightSizeGroup = NULL, intercept = TRUE, loss = c("ls", "logit"), verbose = FALSE, ...) {
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")

  y <- model.response(mf, "numeric")
  X <- model.matrix(mt, mf)
  X <- as.matrix(X)

  res <- MLGL(X, y, hc, lambda, weightLevel, weightSizeGroup, intercept, loss, verbose, ...)

  return(res)
}


#' Hierarchical Clustering with distance matrix computed using bootstrap replicates
#'
#' @param X data
#' @param frac fraction of sample used at each replicate
#' @param B number of replicates
#' @param method desired method: "single", "complete", "average", "mcquitty", "ward.D", "ward.D2", "centroid", "median".
#' @param nCore number of cores
#'
#' @return An object of class \code{hclust}
#'
#'
#' @examples
#' hc <- bootstrapHclust(USArrests, nCore = 1)
#' @export
bootstrapHclust <- function(X, frac = 1, B = 50, method = "ward.D2", nCore = NULL) {
  t1 <- proc.time()
  n <- nrow(X)

  if (frac <= 0 | frac > 1) {
    stop("frac must be between 0 and 1.")
  }


  nInd <- floor(n * frac)
  d <- 0
  for (i in 1:B)
  {
    ind <- sample(n, nInd, replace = TRUE)
    d <- d + parDist(t(X[ind, ]), threads = nCore)
  }

  d <- d / B

  hc <- fastcluster::hclust(d, method = ifelse(is.character(method), method, "ward.D2"))
  t2 <- proc.time()
  tcah <- t2 - t1

  hc$tcah <- as.numeric((t2 - t1)[3])

  return(hc)
}


#
# compute the mimimum weight of each group
#
# @param hc outup of hclust function
#
levelMinWeight <- function(hc, weightLevel = NULL) {
  p <- length(hc$order)

  # highest level at which cluster are seen for the last time
  # the p first are the single variables and the p-2 next are the cluster in the order of apparition
  lvSingle <- sapply((-1):(-p), FUN = function(i) {
    which(hc$merge == i) %% (p - 1)
  })
  lvCluster <- sapply(1:(p - 2), FUN = function(i) {
    which(hc$merge == i) %% (p - 1)
  })
  lvCluster[lvCluster == 0] <- p - 1
  lvCluster <- c(lvCluster, p)

  # branch length. The first one is associated with the partition in 2 clusters
  if (is.null(weightLevel)) {
    weightLevel <- c(0, sqrt(1 / diff(hc$height)), 0)
  }

  # minimum weight of levels of each cluster
  minLevelWeight <- rep(0, 2 * p - 1)
  # minimal weight of single variable
  minLevelWeight[1:p] <- sapply(lvSingle, FUN = function(i) {
    ind <- (weightLevel[1:i] != 0)
    ifelse(sum(ind), min(weightLevel[1:i][ind]), 0)
  }) # If there is only 0, we return 0, else we return the min > 0
  # mininmal weight for groups of 2 and more variables
  minLevelWeight[(p + 1):(2 * p - 1)] <- sapply(1:length(lvCluster), FUN = function(i) {
    ind <- (weightLevel[(i + 1):lvCluster[i]] != 0)
    ifelse(sum(ind), min(weightLevel[(i + 1):lvCluster[i]][ind]), 0)
  })

  return(minLevelWeight)
}



#' Compute the group size weight vector with an authorized maximal size
#'
#' @param hc outup of hclust
#' @param sizeMax maximum size of cluster to consider
#'
#' @return the weight vector
#'
#' @examples
#' set.seed(42)
#' # Simulate gaussian data with block-diagonal variance matrix containing 12 blocks of size 5
#' X <- simuBlockGaussian(50, 12, 5, 0.7)
#' # Generate a response variable
#' y <- X[, c(2, 7, 12)] %*% c(2, 2, -2) + rnorm(50, 0, 0.5)
#' # use 20 as the maximal number of group
#' hc <- hclust(dist(t(X)))
#' w <- computeGroupSizeWeight(hc, sizeMax = 20)
#' # Apply MLGL method
#' res <- MLGL(X, y, hc = hc, weightSizeGroup = w)
#' @export
computeGroupSizeWeight <- function(hc, sizeMax = NULL) {
  uni <- uniqueGroupHclust(hc)
  weight <- as.vector(table(uni$indexGroup))

  if (!is.null(sizeMax)) {
    weight[weight > sizeMax] <- 0
  }

  weight <- sqrt(weight)

  return(weight)
}


#
# @param hc output of hierarchical clustering
#
# @return A matrix with 2 rows, the first row contains the level at which appears each group during the hierarchical clustering
# the second row contains the last level where the group is present. The p first columns represent single variable, the other the cluster in the order
# they appear in the hierarchical clustering
#
levelGroupHC <- function(hc) {
  # Number of variables in the HC
  p <- nrow(hc$merge) + 1

  # Output matrix
  startend <- matrix(nrow = 2, ncol = 2 * p - 1)
  # Level where first appeared each group (j = level containg j groups)
  startend[1, ] <- c(rep(p, p), (p - 1):1)

  # Find the level where each group diassapear
  for (i in 1:nrow(hc$merge))
  {
    for (j in 1:2)
    {
      # a negative number indicates a single variable
      if (hc$merge[i, j] < 0) {
        startend[2, abs(hc$merge[i, j])] <- p - i + 1
      }
      else # a positive number indicates a cluster of 2 or more variables
      {
        startend[2, p + hc$merge[i, j]] <- p - i + 1
      }
    } # end for col os hs$merge
  } # end for row of hc$merge

  # Last group containing all variables
  startend[2, ncol(startend)] <- 1

  rownames(startend) <- c("start", "end")

  return(startend)
}


#
# preliminary step for MLGL. Compute weight, active variables and groups
#
preliminaryStep <- function(hc, weightLevel = NULL, weightSizeGroup = NULL, sizeGroupMax = NULL) {
  # find unique groups of the hclust output
  uni <- uniqueGroupHclust(hc)

  ######## Compute weights
  # compute the minimal weight of partition
  weightLevelGroup <- levelMinWeight(hc, weightLevel)

  # CORRECTION : If weight is infinite, we change in 0 and it will be ignored
  weightLevelGroup[which(is.infinite(weightLevelGroup))] <- 0

  # weight for group size
  if (is.null(weightSizeGroup)) {
    weightSizeGroup <- as.vector(sqrt(table(uni$indexGroup)))
  }

  if (!is.null(sizeGroupMax)) {
    weightSizeGroup[weightSizeGroup > sqrt(sizeGroupMax)] <- 0
  }

  # weight for each group
  weight <- weightSizeGroup * weightLevelGroup

  # new weight without ignored groups
  weightb <- weight
  ignoredGroup <- which(weight == 0) # groups with 0 weights
  weightb <- weightb[-ignoredGroup] # we delete zeros weight

  # level of hc associated to groups
  p <- length(hc$order)
  lv <- levelGroupHC(hc)
  lv <- lv[, -ignoredGroup]

  ######## Create data for gglasso
  varToDelete <- uni$indexGroup %in% ignoredGroup
  var <- uni$varGroup[!varToDelete]
  group <- uni$indexGroup[!varToDelete]

  # group must be consecutively numbered 1,2,3,...
  # need a correction when some groups have to be ignored
  if (length(ignoredGroup) > 0) {
    difNumber <- rep(0, length(group))
    for (i in 1:length(ignoredGroup))
    {
      ind <- which(group > ignoredGroup[i])
      difNumber[ind] <- difNumber[ind] - 1
    }
    group <- group + difNumber
  }

  return(list(group = group, var = var, weight = weightb, level = lv))
}

# check parameters of MLGL function
.checkParameters <- function(X, y, hc, lambda, weightLevel, weightSizeGroup, intercept, verbose, loss, sizeMaxGroup) {
  # check X
  if (!is.matrix(X)) {
    stop("X has to be a matrix.")
  }
  if (any(is.na(X))) {
    stop("Missing values in X not allowed.")
  }
  if (!is.numeric(X)) {
    stop("X has to be a matrix of real.")
  }

  # check y
  if (!is.numeric(y)) {
    stop("y has to be a vector of real.")
  }
  if (any(is.na(y))) {
    stop("Missing values in y not allowed.")
  }
  if (loss == "logit" && any(y %in% c(-1, 1) == FALSE)) {
    stop("Classification method requires the response y to be in {-1,1}")
  }

  # check if X and y are compatible
  if (nrow(X) != length(drop(y))) {
    stop("The length of y and the number of rows of X don't match.")
  }

  # check hc
  if (!is.null(hc)) {
    if (is.character(hc)) {
      if (!(hc %in% c("single", "complete", "average", "mcquitty", "ward.D", "ward.D2", "centroid", "median"))) {
        stop("In character mode, hc must be \"single\", \"complete\", \"average\", \"mcquitty\", \"ward.D\", \"ward.D2\", \"centroid\" or \"median\".")
      }

      if (!is.null(weightLevel)) {
        stop("weightLevel requires a computed hc")
      }
      if (!is.null(weightSizeGroup)) {
        stop("weightSizeGroup requires a computed hc")
      }
    } else {
      # check if hc is a hclust object
      if (class(hc) != "hclust") {
        stop("hc must be an hclust object.")
      }
      # check if hc and X are compatible
      if (length(hc$order) != ncol(X)) {
        stop("hc is not a clustering of the p covariates of X.")
      }

      if (!is.null(weightLevel) && length(weightLevel) != 2 * ncol(X) - 1) {
        stop("weightLevel must be of size 2*p-1")
      }
      if (!is.null(weightSizeGroup) && length(weightSizeGroup) != 2 * ncol(X) - 1) {
        stop("weightSizeGroup must be of size 2*p-1")
      }
    }
  } else {
    if (!is.null(weightLevel)) {
      stop("weightLevel requires the hc argument")
    }
    if (!is.null(weightSizeGroup)) {
      stop("weightSizeGroup requires the hc argument")
    }
  }

  # check if lambda is a vector of positive real
  if (!is.null(lambda)) {
    if (!is.numeric(lambda)) {
      stop("lambda must be a vector of positive real.")
    }
    if (any(lambda < 0)) {
      stop("lambda must be a vector of positive real.")
    }
  }

  # check if weightLevel is a vector of positive real
  if (!is.null(weightLevel)) {
    if (!is.numeric(weightLevel)) {
      stop("weightLevel must be a vector of positive real.")
    }
    if (length(weightLevel) != ncol(X)) {
      stop("weightLevel must have the same length as the number of columns of matrix X.")
    }
    if (any(weightLevel < 0)) {
      stop("weightLevel must be a vector of positive real.")
    }
  }

  # check if weightSizeGroup is a vector of positive real
  if (!is.null(weightSizeGroup)) {
    if (!is.numeric(weightSizeGroup)) {
      stop("weightSizeGroup must be a vector of real.")
    }
    if (any(weightSizeGroup < 0)) {
      stop("weightSizeGroup must be a vector of positive real.")
    }
  }

  # check if intercept is a boolean
  if (length(intercept) != 1) {
    stop("intercept must be a boolean.")
  }
  if (!is.logical(intercept)) {
    stop("intercept must be a boolean.")
  }

  # check if verbose is a boolean
  if (length(verbose) != 1) {
    stop("verbose must be a boolean.")
  }
  if (!is.logical(verbose)) {
    stop("verbose must be a boolean.")
  }

  # check if sizeMaxGroup is a positive integer
  if (!is.null(sizeMaxGroup)) {
    if (length(sizeMaxGroup) != 1) {
      stop("sizeMaxGroup must be a positive integer.")
    }
    if (!.is.wholenumber(sizeMaxGroup)) {
      stop("sizeMaxGroup must be a positive integer.")
    }
  }


  invisible(return(NULL))
}
