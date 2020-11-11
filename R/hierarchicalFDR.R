
#### intern function
# method from hierarchical false discovery rate-controlling method
# fdr control of famlily : q
# fdr control of full tree : q * delta * 2  (delta < 1.44)
# fdr control of outer node : q * L * delta * 2 (delta < 1.44)
hierarchicalFDRTesting <- function(hierMat, group, grouplm, X, y, test = partialFtest) {
  # root of the tree
  indRoot <- findRoot2(hierMat)

  # output container
  pvalues <- rep(0, nrow(hierMat))
  adjPvalues <- rep(0, nrow(hierMat))

  # all Leaves in the hierarchical matrix
  allLeaves <- which(rowSums(hierMat) == 1)

  familyToTest <- list()
  familyToTest[[1]] <- indRoot
  continue <- TRUE

  # hierarchical testing
  while (continue) {
    familyToTestTemp <- list()
    compteur <- 1
    continue <- FALSE
    # current family at the level
    for (i in 1:length(familyToTest))
    {
      # test the group of a family
      for (gr in familyToTest[[i]])
      {
        # group included in gr
        subGroup <- setdiff(which(hierMat[gr, ]), gr)

        # if no subgroups, gr is a leaf
        if (length(subGroup) == 0) {
          subGroup <- gr
        }

        # only group corresponding to leaves
        subLeaves <- intersect(subGroup, allLeaves)
        toTest0 <- which(grouplm %in% group[subLeaves])
        pvalues[gr] <- test(X, y, toTest0)
      }

      # BH correction for the family
      adjPvalues[familyToTest[[i]]] <- p.adjust(pvalues[familyToTest[[i]]], "BH")

      # if selected, we look for the children to test
      for (j in 1:length(familyToTest[[i]]))
      {
        child <- children(familyToTest[[i]][j], hierMat)
        if (length(child) > 0) {
          familyToTestTemp[[compteur]] <- child
          compteur <- compteur + 1
          continue <- TRUE
        }
      } # end for selected group
    } # end for family
    familyToTest <- familyToTestTemp
  } # end while hierarchy


  return(list(pvalues = pvalues, adjPvalues = adjPvalues, groupId = as.numeric(colnames(hierMat))))
}


#'
#' Apply hierarchical test for each hierarchy, and test external variables for FDR control at level alpha
#'
#' @title Hierachical testing with FDR control
#'
#' @param X original data
#' @param y associated response
#' @param group vector with index of groups. group[i] contains the index of the group of the variable var[i].
#' @param var vector whith the variables contained in each group. group[i] contains the index of the group of the variable var[i].
#' @param test function for testing the nullity of a group of coefficients in linear regression. 3 parameters : X : design matrix, y response and varToTest : vector of variables to test; return a pvalue
#'
#' @return a list containing :
#' \describe{
#'   \item{pvalues}{pvalues of the different test (without correction)}
#'   \item{adjPvalues}{adjusted pvalues}
#'   \item{groupId}{Index of the group}
#'   \item{hierMatrix}{Matrix describing the hierarchical tree.}
#'   }
#'
#' @details
#' Version of the hierarchical testing procedure of Yekutieli for MLGL output. You can use th \link{selFDR} function to select groups
#' at a desired level alpha.
#'
#'
#' @examples
#' set.seed(42)
#' X <- simuBlockGaussian(50, 12, 5, 0.7)
#' y <- X[, c(2, 7, 12)] %*% c(2, 2, -2) + rnorm(50, 0, 0.5)
#' res <- MLGL(X, y)
#' test <- hierarchicalFDR(X, y, res$group[[20]], res$var[[20]])
#' @references Yekutieli, Daniel. "Hierarchical False Discovery Rate-Controlling Methodology." Journal of the American Statistical Association 103.481 (2008): 309-16.
#'
#' @seealso \link{selFDR}, \link{hierarchicalFWER}
#'
#' @export
hierarchicalFDR <- function(X, y, group, var, test = partialFtest) {
  # check parameters
  .checkhierarchicalFWER(X, y, group, var, test, TRUE)

  # hierarchical matrix
  hierInfo <- groupHier(group, var)

  # lm with leaves represented by their first principal component
  reslm <- acpOLStest(X, y, hierInfo$grouplm, hierInfo$varlm)

  # new hierMat with complementary group
  #   hierMatTot <- compHierMatTot(hierInfo)
  hierMatTot <- hierInfo$hierTot
  groupId <- colnames(hierMatTot)

  # hierarchical testing on the tree of indRoot
  out <- hierarchicalFDRTesting(hierMatTot, groupId, reslm$group, reslm$newdata, y, test = partialFtest)

  out$group <- hierInfo$groupTot
  out$var <- hierInfo$varTot
  out$hierMatrix <- hierMatTot

  return(out)
}


#'
#' Select groups from hierarchical testing procedure with FDR control (\link{hierarchicalFDR})
#'
#' @title Selection from hierarchical testing with FDR control
#'
#' @param out output of \link{hierarchicalFDR} function
#' @param alpha control level for test
#' @param global if FALSE the provided alpha is the desired level control for each family.
#' @param outer if TRUE, the FDR is controlled only on outer node (rejected groups without rejected children) . If FALSE, it is controlled on the full tree.
#'
#' @return a list containing :
#' \describe{
#'   \item{toSel}{vector of boolean. TRUE if the group is selected}
#'   \item{groupId}{Names of groups}
#'   \item{local.alpha}{control level for each family of hypothesis}
#'   \item{global.alpha}{control level for the tree (full tree or outer node)}
#'   }
#'
#' @details
#' See the reference for mode details about the method.
#'
#' If each family is controlled at a level alpha, we have the following control :
#' FDR control of full tree : alpha * delta * 2  (delta = 1.44)
#' FDR control of outer node : alpha * L * delta * 2 (delta = 1.44)
#'
#' @examples
#' set.seed(42)
#' X <- simuBlockGaussian(50, 12, 5, 0.7)
#' y <- X[, c(2, 7, 12)] %*% c(2, 2, -2) + rnorm(50, 0, 0.5)
#' res <- MLGL(X, y)
#' test <- hierarchicalFDR(X, y, res$group[[20]], res$var[[20]])
#' sel <- selFDR(test, alpha = 0.05)
#' @references Yekutieli, Daniel. "Hierarchical False Discovery Rate-Controlling Methodology." Journal of the American Statistical Association 103.481 (2008): 309-16.
#'
#' @seealso \link{hierarchicalFDR}
#'
#' @export
selFDR <- function(out, alpha = 0.05, global = TRUE, outer = TRUE) {
  ## check arguments
  .checkselFDR(out$adjPvalues, out$hierMatrix, alpha, global, outer)

  delta <- 1.44
  # number of levels of the hierarchy
  L <- ifelse(outer, numberLevels(out$hierMatrix), 1)
  # compute local and global alpha
  local.alpha <- ifelse(global, alpha / (delta * 2 * L), alpha)
  global.alpha <- ifelse(global, alpha, min(alpha * (delta * 2 * L), 1))
  # if global = 1, then use local = 1
  local.alpha <- ifelse(global.alpha == 1, 1, local.alpha)
  # if one level, simphe BH adjust
  local.alpha <- ifelse(outer & (L == 1), alpha, local.alpha)
  global.alpha <- ifelse(outer & (L == 1), alpha, global.alpha)

  #
  toSel <- rep(FALSE, length(out$adjPvalues))

  # indice of groupes at the top of hierarchy
  family <- findRoot2(out$hierMatrix)

  continue <- TRUE

  while (continue) {
    continue <- FALSE
    # select group with adjusted pavalues <= local.alpha
    toSel[family] <- (out$adjPvalues[family] <= local.alpha)

    # find children of selected group
    if (any(toSel[family])) {
      ind <- which(toSel[family])
      newfamily <- c()
      for (i in 1:length(ind)) {
        newfamily <- c(newfamily, children(family[ind[i]], out$hierMatrix))
      }

      family <- newfamily
      continue <- (length(family) > 0)
    }
  } # fin while

  # we select only outer node
  if (outer) {
    toSel[which(toSel)[which(rowSums(out$hierMatrix[toSel, toSel, drop = FALSE]) > 1)]] <- FALSE
  }

  return(list(toSel = toSel, groupId = as.numeric(colnames(out$hierMatrix)), local.alpha = local.alpha, global.alpha = global.alpha))
}


## check paraleters of selFWER function
.checkselFDR <- function(adjPvalues, hierMatrix, alpha, global, outer) {
  .checkselFWER(adjPvalues, hierMatrix, alpha)

  # check if global is a boolean
  if (length(global) != 1) {
    stop("global must be a boolean.")
  }
  if (!is.logical(global)) {
    stop("global must be a boolean.")
  }

  # check if outer is a boolean
  if (length(outer) != 1) {
    stop("outer must be a boolean.")
  }
  if (!is.logical(outer)) {
    stop("outer must be a boolean.")
  }

  invisible(return(NULL))
}
