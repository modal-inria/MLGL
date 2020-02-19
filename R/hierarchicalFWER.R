
#
# Computation of the Shaffer coefficient
#
# pvalue correction if all siblings of a group are leaves nodes
#
# ind : index of the tested group
# hierMatTot : matrix describing the hierarchy
# return sh : shaffer coefficient (apply before hierarchical adjustment)
shafferImprovement <- function(ind, hierMatTot)
{
  # size of group in (in number of leaves)
  sizeInd = sum(hierMatTot[ind,which(rowSums(hierMatTot)==1)])
  # hierarchy level of group ind : 1 = top
  levelHier = colSums(hierMatTot)
  
  # ancestor of group ind
  parents = setdiff(which(hierMatTot[,ind]),ind)
  # parent of group in
  parent = parents[which.min(colSums(hierMatTot[,parents,drop=FALSE]))]
  # siblings of group ind
  sibling = intersect(setdiff(which(hierMatTot[parent,]),ind),which(levelHier==levelHier[ind]))
  nSibling = length(sibling)
  # Number of sibling with descendants
  nSiblingdiff = length(which(rowSums(hierMatTot[sibling,,drop = FALSE]) > 1))
  
  sh = 1
  # We can't apply the Shaffer correction if all siblings have not any descendants
  if( (nSiblingdiff == 0) & (nSibling > 0) )
    sh = (sizeInd + nSibling)/(sizeInd + nSibling -1)
  
  return(sh)
}

#
# hierarchical ajustment for p-values in the hierarchical testing procedure
#
# pvalue : pvalue of the test 
# gr : index of group in the hierchical matrix
# hierMat : hierrarchical matrix
# adjPvalues : all pvalues already computed
# Shaffer : boolean, if TRUE, a shaffer correction is performed
hierarchicalAdjustement <- function(pvalue, gr, hierMat, adjPvalues, sizeRoot, sizeGr, Shaffer = FALSE)
{
  # Shaffer adjustement
  sh <- ifelse(Shaffer, shafferImprovement(gr, hierMat), 1)
  
  adjPvalue = min(pvalue * sizeRoot / sizeGr / sh, 1) # Meinshausen's group size adjustement
  adjPvalue = max(c(adjPvalues[which(hierMat[,gr])], adjPvalue), na.rm = TRUE) # Meinshausen's hierarchical adjustement
  
  return(adjPvalue)  
}

# hierarchical testing for a completed tree 
#
# indRoot ; index of the root of the tree in the hierarchy matrix
# hierMat : matrix describing the hierarchy
# group : name of groups
# grouplm : group used in OLS for the test function (leaves of the tree)
# X : matrix with principal component
# y : response
# test : function used for test
# Shaffer : apply the Shaffer correction ?
#
#
hierarchicalTesting <- function(indRoot, hierMat, group, grouplm, X, y, test = partialFtest, Shaffer = FALSE)
{
  # all Leaves in the hierarchical matrix
  a <- (rowSums(hierMat) == 1)
  allLeaves <- which(a)
  indcont   <- which(!a)
  
  # all group included in the root indRoot
  indGrInHier <- which(hierMat[indRoot,])
  grInHier <- group[indGrInHier]
  
  # Only the leaves of the hierarchy with root indRoot    
  grRoot   <- intersect(indGrInHier, allLeaves)
  sizeRoot <- length(grRoot) # number of leaves in the tree
  
  # Number of groups (leaves and no leaves in the hierarchy)
  sizeHier <- length(grInHier)
  
  # output
  pvalues = adjPvalues <- rep(NA, sizeHier)
  
  continue = TRUE
  
  # Index of the non rejected group to test at the current level
  indToTest <- indRoot
  
  #hierarchical testing
  while(continue)
  {  
    continue = FALSE
    indToTestTemp = c()
    
    #for each  group of the current level
    for(gr in indToTest)
    {
      # group included in gr
      subGroup <- group[setdiff(which(hierMat[gr,]),gr)]
      if(length(subGroup)==0)
        subGroup = group[gr]
      
      # only group corresponding to leaves    
      subLeaves <- intersect(subGroup, group[grRoot])
      sizeGr   <- length(subLeaves)
      toTest0 <- which(grouplm %in% subLeaves)
      
      indGrOutObj <- match(group[gr], grInHier)
      pvalues[indGrOutObj]    = test(X, y, toTest0)
      if(length(indGrOutObj)>1)
      {
        adjPvalues[indGrOutObj] = hierarchicalAdjustement(pvalues[indGrOutObj], indGrOutObj, hierMat[indGrInHier, indGrInHier], adjPvalues, sizeRoot, sizeGr, Shaffer)
      }
      else
      {
        adjPvalues[indGrOutObj] = pvalues[indGrOutObj]
      }   
    }#end for group
    
    for(j in 1:length(indToTest))
    {
      child = children(indToTest[j], hierMat)
      if(length(child) > 0)
      {
        indToTestTemp = c(indToTestTemp, child)
        continue = TRUE
      }
      
    }# end for selected group
    indToTest = indToTestTemp
    
  }#end while hierarchy
  
  
  return(list(pvalues = pvalues, adjPvalues = adjPvalues))
}


# hierarchicalTestingNew <- function(indRoot, hierMat, group, grouplm, X, y, test = partialFtest, Shaffer = FALSE)
# {
#   # all Leaves in the hierarchical matrix
#   a <- (rowSums(hierMat) == 1)
#   allLeaves <- which(a)
#   indcont   <- which(!a)
#   
#   # all group included in the root indRoot
#   indGrInHier <- which(hierMat[indRoot,])
#   grInHier <- group[indGrInHier]
#   
#   # Only the leaves of the hierarchy with root indRoot    
#   grRoot   <- intersect(indGrInHier, allLeaves)
#   sizeRoot <- length(grRoot) # number of leaves in the tree
#   
#   # Number of groups (leaves and no leaves in the hierarchy)
#   sizeHier <- length(grInHier)
#   
#   # output
#   pvalues = adjPvalues <- rep(NA, sizeHier)
#   
#   continue = TRUE
#   
#   # Index of the non rejected group to test at the current level
#   indToTest <- indRoot
#   
#   
#   sizeGrInLeaves <- rep(0,nrow(grInHier))
#   i <- 1
#   for(gr in indGrInHier)
#   {
#     
#     ### Find the leaves of the group
#     # group included in gr
#     subGroup <- group[setdiff(which(hierMat[gr,]),gr)]
#     if(length(subGroup)==0)
#       subGroup = group[gr]
#     
#     # only group corresponding to leaves    
#     subLeaves <- intersect(subGroup, group[grRoot])
#     
#     sizeGrInLeaves[i] = length(subLeaves) 
#     sizeGr <- sizeGrInLeaves[i]
#     
#     toTest0 <- which(grouplm %in% subLeaves)
#     
#     pvalues[indGrOutObj]    = test(X, y, toTest0)
#     if(length(indGrOutObj)>1)
#     {
#       adjPvalues[indGrOutObj] = hierarchicalAdjustement(pvalues[indGrOutObj], indGrOutObj, hierMat[indGrInHier, indGrInHier], adjPvalues, sizeRoot, sizeGr, Shaffer)
#     }
#     else
#     {
#       adjPvalues[indGrOutObj] = pvalues[indGrOutObj]
#     }   
#   
#     # TODO
#     # ajustement par la taille : ok
#     # ajustement hierarchique a voir : faire une liste ancestor pvaladj[gr] = max(pval[gr],pval[ancestor[[gr]]]) : pval contient ceux deja ajuste par la taille
#     # mettre shaffer = FALSE
#     
#     i <- i + 1
#   }# fin for group hier
#   
# #   
# #   #hierarchical testing
# #   while(continue)
# #   {  
# #     continue = FALSE
# #     indToTestTemp = c()
# #     
# #     #for each  group of the current level
# #     for(gr in indToTest)
# #     {
# #       # group included in gr
# #       subGroup <- group[setdiff(which(hierMat[gr,]),gr)]
# #       if(length(subGroup)==0)
# #         subGroup = group[gr]
# #       
# #       # only group corresponding to leaves    
# #       subLeaves <- intersect(subGroup, group[grRoot])
# #       sizeGr   <- length(subLeaves)
# #       toTest0 <- which(grouplm %in% subLeaves)
# #       
# #       indGrOutObj <- match(group[gr], grInHier)
# #       
# #       pvalues[indGrOutObj]    = test(X, y, toTest0)
# #       if(length(indGrOutObj)>1)
# #       {
# #         adjPvalues[indGrOutObj] = hierarchicalAdjustement(pvalues[indGrOutObj], indGrOutObj, hierMat[indGrInHier, indGrInHier], adjPvalues, sizeRoot, sizeGr, Shaffer)
# #       }
# #       else
# #       {
# #         adjPvalues[indGrOutObj] = pvalues[indGrOutObj]
# #       }   
# #     }#end for group
# #     
# #     for(j in 1:length(indToTest))
# #     {
# #       child = children(indToTest[j], hierMat)
# #       if(length(child) > 0)
# #       {
# #         indToTestTemp = c(indToTestTemp, child)
# #         continue = TRUE
# #       }
# #       
# #     }# end for selected group
# #     indToTest = indToTestTemp
# #     
# #   }#end while hierarchy
# #   
#   
#   return(list(pvalues = pvalues, adjPvalues = adjPvalues))
# }



#'
#' Apply hierarchical test for each hierarchy, and test external variables for FWER control at level alpha   
#'
#' @title Hierachical testing with FWER control 
#'
#' @param X original data
#' @param y associated response
#' @param group vector with index of groups. group[i] contains the index of the group of the variable var[i].
#' @param var vector whith the variables contained in each group. group[i] contains the index of the group of the variable var[i].
#' @param test function for testing the nullity of a group of coefficients in linear regression. 3 parameters : X : design matrix, y response and varToTest : vector of variables to test; return a pvalue
#' @param Shaffer boolean, if TRUE, a shaffer correction is performed
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
#' Version of the hierarchical testing procedure of Meinshausen for MLGL output. You can use th \link{selFWER} function to select groups
#' at a desired level alph
#' 
#' 
#' @examples 
#' set.seed(42)
#' X = simuBlockGaussian(50,12,5,0.7)
#' y = X[,c(2,7,12)]%*%c(2,2,-2) + rnorm(50,0,0.5)
#' res = MLGL(X,y)
#' test = hierarchicalFWER(X, y, res$group[[20]], res$var[[20]])
#' 
#' 
#' @references Meinshausen, Nicolai. "Hierarchical Testing of Variable Importance." Biometrika 95.2 (2008): 265-78.
#' 
#' @seealso \link{selFWER}, \link{hierarchicalFDR}
#' 
#' @export
hierarchicalFWER <- function(X, y, group, var, test = partialFtest, Shaffer = FALSE)
{
  # check parameters
  .checkhierarchicalFWER(X, y, group, var, test, Shaffer)
    
  
  # hierarchical matrix
  hierInfo <- groupHier(group, var)
  grdif    <- unique(hierInfo$groupTot)
  grouplm  <- unique(hierInfo$grouplm)
  # total number of leaves
  m <- length(grouplm)
  # lm with leaves represented by their first principal component
  reslm <- acpOLStest(X, y, hierInfo$grouplm, hierInfo$varlm)

  # new hierMat with complementary group
  hierMatTot <- hierInfo$hierTot    
  
  # indice of groups at the top of a hierarchy 
  indGrTop <- findRoot2(hierMatTot)
  
  # output container
  pvalues    <- rep(NA, length(grdif))
  adjPvalues <- rep(NA, length(grdif))
  
  ## Hierarchical test on the different trees
  for(indRoot in indGrTop)
  {
    # hierarchical testing on the tree of indRoot
    out <- hierarchicalTesting(indRoot, hierMatTot, grdif, grouplm, reslm$newdata, reslm$y, test, Shaffer)
    
    # index of group contained in the hierarchy
    indGrHierRoot <- which(hierMatTot[indRoot,])
    
    # indice of leaves and single variables
    indLeaves <- leaves(hierMatTot[indGrHierRoot, indGrHierRoot, drop = FALSE])
    
    # stock results in output
    pvalues[indGrHierRoot]    = out$pvalues
    adjPvalues[indGrHierRoot] = pmin(out$adjPvalues * m/length(indLeaves), 1)# adjustement for multiple tree # each tree is penalized by its size
    
  }# end for root
  
  groupId <- as.numeric(rownames(hierMatTot)) 
  return(list(pvalues = pvalues, adjPvalues = adjPvalues,  groupId = groupId, hierMatrix = hierMatTot))
}

#'
#' Select groups from hierarchical testing procedure with FWER control (\link{hierarchicalFWER})
#'
#' @title Selection from hierarchical testing with FWER control
#'
#' @param out output of \link{hierarchicalFDR} function
#' @param alpha control level for test
#' 
#' @return a list containing :
#' \describe{
#'   \item{toSel}{vector of boolean. TRUE if the group is selected}
#'   \item{groupId}{Names of groups}
#'   }
#'
#' @details 
#' Only outer nodes (rejected groups without rejected children) are returned as TRUE.
#' 
#' 
#' @examples 
#' set.seed(42)
#' X = simuBlockGaussian(50,12,5,0.7)
#' y = X[,c(2,7,12)]%*%c(2,2,-2) + rnorm(50,0,0.5)
#' res = MLGL(X,y)
#' test = hierarchicalFWER(X, y, res$group[[20]], res$var[[20]])
#' sel = selFWER (test, alpha = 0.05)
#' 
#' 
#' @references Meinshausen, Nicolai. "Hierarchical Testing of Variable Importance." Biometrika 95.2 (2008): 265-78.
#' 
#' @seealso \link{hierarchicalFWER}
#' 
#' @export
selFWER <- function(out, alpha = 0.05)
{
  # check parameters
  .checkselFWER(out$adjPvalues, out$hierMatrix, alpha)
  
  # output
  toSel = rep(FALSE, length(out$adjPvalues))
  
  # indice of groupes at the top of hierarchy
  family = findRoot2(out$hierMatrix)
  
  continue = TRUE
  
  while(continue)
  {
    continue = FALSE
    # select group with adjusted pavalues <= local.alpha
    toSel[family] = (out$adjPvalues[family] <= alpha)
    
    # find children of selected group
    if(any(toSel[family]))
    {
      ind = which(toSel[family])
      newfamily = c()
      for(i in 1:length(ind))
        newfamily = c(newfamily, children(family[ind[i]], out$hierMatrix))
      
      family = newfamily
      continue = (length(family)>0)
    } 
    
  }#fin while
  
  # we select only outer node
  toSel[which(toSel)[which(rowSums(out$hierMatrix[toSel,toSel,drop=FALSE])>1)]] = FALSE
  
  
  return(list(toSel = toSel, groupId = as.numeric(colnames(out$hierMatrix))))
}

## check parameters of hierarchicalFWER function
.checkhierarchicalFWER <- function(X, y, group, var, test, Shaffer)
{
  
  #check X
  if(!is.matrix(X)) 
    stop("X has to be a matrix.")
  if(any(is.na(X))) 
    stop("Missing values in X not allowed.")
  if(!is.numeric(X))
    stop("X has to be a matrix of real.")
  
  #check y
  if(!is.numeric(y))
    stop("y has to be a vector of real.")
  if(any(is.na(y))) 
    stop("Missing values in y not allowed.")
  
  #check if X and y are compatible
  if(nrow(X)!=length(drop(y)))
    stop("The length of y and the number of rows of X don't match.")
  
  #check group
  if(!is.numeric(group))
    stop("group has to be a vector of integer.")
  if(any(is.na(group))) 
    stop("Missing values in group not allowed.")  

  #check var
  if(!is.numeric(var))
    stop("var has to be a vector of integer.")
  if(any(is.na(var))) 
    stop("Missing values in var not allowed.")
  if(any(var>ncol(X)))
    stop("var must contains index of column of X.")
  if(any(var <= 0))
    stop("var must be a vector of positive integer.")
  
  #check if group and var are compatible
  if(length(drop(group))!=length(drop(var)))
    stop("The length of group and var don't match.")
  
  # check if test is a function
  if(!is.function(test))
    stop("test must be a function.")
  
  #check if Shaffer is a boolean
  if(length(Shaffer)!=1)
    stop("Shaffer must be a boolean.")
  if(!is.logical(Shaffer))
    stop("Shaffer must be a boolean.")
  
  invisible(return(NULL))
}


## check parameters of selFWER function
.checkselFWER <- function(adjPvalues, hierMatrix, alpha)
{
  #alpha
  if(length(alpha)!=1)
    stop("alpha must be a real between 0 and 1.")
  if((alpha <=0) || (alpha>1))
    stop("alpha must be a real between 0 and 1.")
 
  # hierMatrix
  if(!is.matrix(hierMatrix)) 
    stop("hierMatrix has to be a matrix.")
  if(any(is.na(hierMatrix))) 
    stop("Missing values in hierMatrix not allowed.")
  if(!is.logical(hierMatrix))
    stop("hierMatrix has to be a matrix of boolean.")

  # adjPvalues
  if(!is.numeric(adjPvalues))
    stop("adjPvalues has to be a vector of real between 0 and 1.")

  # if(any(is.na(adjPvalues))) 
  #   stop("Missing values in adjPvalues not allowed.")
  if(any(adjPvalues[!is.na(adjPvalues)] <= 0) || any(adjPvalues[!is.na(adjPvalues)]  > 1))
    stop("adjPvalues has to be a vector of real between 0 and 1.")

  invisible(return(NULL))
}
