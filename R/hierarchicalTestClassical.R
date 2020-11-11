
hierarchicalTestingClassical <- function(X, y, hc = NULL, alpha = 0.05, test = partialFtest) {
  # perform hierarchical clustering if not provided
  if (is.null(hc)) {
    Xb <- scale(X)
    d <- dist(t(Xb))
    hc <- fastcluster::hclust(d, method = "ward.D2")
  }

  # find unique group
  uniqueGroup <- uniqueGroupHclust(hc)


  group <- uniqueGroup$indexGroup
  var <- uniqueGroup$varGroup

  # information about hierarchy
  hierInfo <- groupHier(group, var)

  grdif <- unique(hierInfo$groupTot)
  # leaves of the tree
  grouplm <- unique(hierInfo$grouplm)

  # matrix describieng the tree
  hierMatTot <- hierInfo$hierTot

  # indice of groups at the top of a hierarchy
  indGrTop <- findRoot2(hierMatTot)

  # hierarchical testing
  out <- hierarchicalTesting(indGrTop, hierMatTot, grdif, grouplm, X, y, partialFtest, TRUE)
  out$hierMatrix <- hierMatTot

  # return selected groups
  outSel <- selFWER(out, alpha = alpha)

  return(c(out, outSel, list(group = group, var = var, alpha = alpha, hierMat = hierMatTot)))
}



hierarchicalFDRTestingClassical <- function(X, y, hc = NULL, alpha = 0.05, test = partialFtest, global = TRUE, outer = TRUE) {
  # perform hierarchical clustering if not provided
  if (is.null(hc)) {
    Xb <- scale(X)
    d <- dist(t(Xb))
    hc <- fastcluster::hclust(d, method = "ward.D2")
  }

  # find unique group
  uniqueGroup <- uniqueGroupHclust(hc)


  group <- uniqueGroup$indexGroup
  var <- uniqueGroup$varGroup

  # information about hierarchy
  hierInfo <- groupHier(group, var)

  grdif <- unique(hierInfo$groupTot)
  # leaves of the tree
  grouplm <- unique(hierInfo$grouplm)

  # matrix describieng the tree
  hierMatTot <- hierInfo$hierTot

  # indice of groups at the top of a hierarchy
  indGrTop <- findRoot2(hierMatTot)

  # hierarchical testing
  out <- hierarchicalFDRTesting(hierMatTot, grdif, grouplm, X, y, test = test)
  out$hierMatrix <- hierMatTot

  # return selected groups
  outSel <- selFDR(out, alpha = alpha, global = global, outer = outer)

  return(c(out, outSel, list(group = group, var = var, hierMat = hierMatTot)))
}
