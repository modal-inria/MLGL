# find root of the tree
findRoot <- function(hierMatrix) {
  # indice of groups containing other groups
  indcont <- which(rowSums(hierMatrix) != 1)
  # indice of groups not included in a group
  indtop <- which(colSums(hierMatrix) == 1)

  # indice of groups at the top of a hierarchy
  indGrTop <- intersect(indtop, indcont)

  return(indGrTop)
}

# find roots of trees and isolated groups
findRoot2 <- function(hierMatrix) {
  # indice of groups not included in a group
  indtop <- which(colSums(hierMatrix) == 1)
  names(indtop) <- NULL
  return(indtop)
}

# find isolated groups (not belongging to hierarchical trees)
findIsolateVariable <- function(hierMatrix) {
  # indice of groups containing other groups
  indcont <- which(rowSums(hierMatrix) == 1)
  # indice of groups not included in a group
  indtop <- which(colSums(hierMatrix) == 1)

  # indice of groups at the top of a hierarchy
  indGrTop <- intersect(indtop, indcont)

  return(indGrTop)
}

# find direct subgroups of a group
children <- function(ind, hierMatrix) {
  # all descendants
  descendant <- setdiff(which(hierMatrix[ind, ]), ind)

  children <- descendant
  if (length(descendant) > 1) {
    # we delete descendant of descendant
    for (i in descendant)
    {
      children <- setdiff(children, setdiff(which(hierMatrix[i, ]), i))
    }
  }

  return(children)
}

# find parent of a group
parent <- function(ind, hierMatrix) {
  # all ancestors
  ancestor <- setdiff(which(hierMatrix[, ind]), ind)

  parent <- ancestor
  if (length(ancestor) > 1) {
    # we delete ancestor of ancestor
    for (i in ancestor)
    {
      parent <- setdiff(parent, setdiff(which(hierMatrix[, i]), i))
    }
  }

  return(parent)
}

# retrun index of the leaves of the tree (and single variables)
leaves <- function(hierMatrix) {
  return(which(rowSums(hierMatrix[, , drop = FALSE]) == 1))
}

# compute the depth of a tree : the maximum number of levels from root to leaves
numberLevels <- function(hierMatrix) {

  # leaves of the tree
  allLeaves <- leaves(hierMatrix)

  L <- rep(0, length(allLeaves))

  for (i in seq_along(L))
  {
    ok <- TRUE
    ind <- allLeaves[i]
    while (ok) {
      ind <- parent(ind, hierMatrix)
      ok <- (length(ind) > 0)
      L[i] <- L[i] + 1
    }
  }

  return(max(L))
}

# get only selected outer nodes (groups do not containing rejected children) from a selection
outerNode <- function(toSel, hierMatrix) {
  return(rowSums(hierMatrix[toSel, toSel, drop = FALSE]) <= 1)
}

# compute matrix describing the hierarchy
#
# return the hierarchy matrix : a binary square matrix. Each row and each column represents a different group.
# row i col j = TRUE if jth group is included in ith group, FALSE otherwise
#
#
compHierMat <- function(group, var) {
  grdif <- unique(group)
  hierMat <- matrix(FALSE, nrow = length(grdif), ncol = length(grdif))
  j <- 1

  grdif <- sort(grdif)
  for (indgr in grdif)
  {
    indg <- which(group == indgr)
    hierMat[j, ] <- tapply(var, group, FUN = function(x) {
      all(x %in% var[indg])
    })
    j <- j + 1
  }

  colnames(hierMat) <- rownames(hierMat) <- grdif

  return(hierMat)
}

# compute the matrix describing the completed hierarchy
compHierMatTot <- function(hierInfo) {
  varTot <- c(hierInfo$var, hierInfo$varComp)
  groupTot <- c(hierInfo$group, hierInfo$groupComp)

  hierMatTot <- compHierMat(groupTot, varTot)

  return(hierMatTot)
}



# This function gives the group to keep for the hierarchical test procedure
#
# group: vector containing the index of group selected
# var: vector containing the index of variables containing in the different selected groups
# addRoot: if TRUE, add a group containing all the groups
# group and var have the same size
#
# return the hierarchy matrix, the completed hieraechy matrix and vectors describing groups
#
groupHier <- function(group, var, addRoot = FALSE) {
  # unique group selected
  grdif <- unique(group)

  # matrix describing the relation of inclusion between groups
  # (i,j) is TRUE if group grdif[j] is included in grdif[i]
  hierMat <- compHierMat(group, var)


  # groupes contenant d'autres groupes
  rSumHM <- rowSums(hierMat)
  grContGr <- which(rSumHM != 1)

  # is there a root (a group containing all the other group)
  isRoot <- (max(rSumHM) == length(grdif))

  # groupe ne contenant pas d'autres groupes
  grNotCont <- grdif
  if (length(grContGr) != 0) {
    grNotCont <- grdif[-grContGr]
  }


  indToKeep <- group %in% grNotCont
  grouplm <- group[indToKeep]
  varlm <- var[indToKeep]

  # complementary groups : group to add to have a full hierarchy
  GRCOMP <- list()
  groupComp <- varComp <- c()
  if (length(grContGr) > 0) {
    for (i in grContGr)
    {
      grIncGrI <- setdiff(which(hierMat[i, ]), i)
      indGrIncGrI <- which(group %in% grdif[grIncGrI])
      indGrI <- which(group %in% grdif[i])
      varGrComp <- setdiff(var[indGrI], var[indGrIncGrI])

      GRCOMP[[i]] <- varGrComp
    }

    # we add complementary groups
    lcomp <- sapply(GRCOMP, length)

    if (any(lcomp > 0)) {
      groupComp <- -rep(1:length(lcomp[lcomp != 0]), lcomp[lcomp != 0])
      varComp <- unlist(GRCOMP)
    }
  }

  grouplm <- c(grouplm, groupComp)
  varlm <- c(varlm, varComp)

  groupTot <- c(group, groupComp)
  varTot <- c(var, varComp)

  # root
  if (!isRoot & addRoot) {
    # unique variables
    varUni <- unique(varTot)

    # add the group to the current list
    numGrRoot <- min(groupComp, 0) - 1
    groupTot <- c(groupTot, rep(numGrRoot, length(varUni)))
    groupComp <- c(groupComp, rep(numGrRoot, length(varUni)))

    # add all the variables to the current list
    varTot <- c(varTot, varUni)
    varComp <- c(varComp, varUni)
  }

  # order group number (need for grouphier and other functions)
  ordGroupTot <- order(groupTot)
  groupTot <- groupTot[ordGroupTot]
  varTot <- varTot[ordGroupTot]

  ordGroupLm <- order(grouplm)
  grouplm <- grouplm[ordGroupLm]
  varlm <- varlm[ordGroupLm]


  if (setequal(groupTot, group)) {
    hierMatTot <- hierMat
  } else {
    hierMatTot <- compHierMat(groupTot, varTot)
  }

  return(list(
    hier = hierMat, group = group, var = var, groupComp = groupComp, varComp = varComp,
    varlm = varlm, grouplm = grouplm, varTot = varTot, groupTot = groupTot, hierTot = hierMatTot
  ))
}
