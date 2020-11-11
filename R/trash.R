# thmSolPath <- function(res, X, y)
# {
#   acp <- list()
#   acpdone <- c()
#
#   REJECT <- list()
#   nbReject <- rep(0, length(res$lambda))
#   prevSelGroup = selGroup <- c()
#   for(i in 1:length(res$lambda))
#   {
#     # if no groups are selected we do nothing
#     if(length(res$group[[i]])>0)
#     {
#       selGroup = unique(res$group[[i]])
#
#       # if the selected groups have not changed compared with the last iteration, we copy the result
#       if(setequal(prevSelGroup, selGroup))
#       {
#         REJECT[[i]] = REJECT[[i-1]]
#         nbReject[i] = nbReject[i-1]
#       }
#       else
#       {
#
#
#
#
#
#
#
#
#         # compute acp
#         indGr <- unique(hierInfo$grouplm)
#
#         indToComp = which(!(indGr%in%done))
#         for(ind in indToComp)
#         {
#           hierInfo$varlm[group%in%ind]
#           acp[[ind]] = PCA(X[,indvar], scale.unit=TRUE, ncp=1)$ind$coord
#         }
#
#         newdata = do.call(cbind,acp[indGr])
#
#
#
#
#
#
#
#
#
#       }# end if no selection
#
#       prevSelGroup = selGroup
#
#     }# end for lambda
#
#
#   }# fin function
