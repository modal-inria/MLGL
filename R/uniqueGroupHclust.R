#' Find all unique groups in \code{\link{hclust}} results
#'
#' @author Quentin Grimonprez
#' @param hc output of \code{\link{hclust}} function
#' @return A list containing:
#' \describe{
#'   \item{indexGroup}{Vector containing the index of variables.}
#'   \item{varGroup}{Vector containing the index of the group of each variable.}
#' }
#'
#' @examples
#' hc <- hclust(dist(USArrests), "average")
#' res <- uniqueGroupHclust(hc)
#' @export
uniqueGroupHclust <- function(hc) {
  # check if hc is a hclust object
  if (class(hc) != "hclust") {
    stop("hc must be an hclust object.")
  }

  # hc$merge contains the order and composition (2 number) of the different merge
  nr <- nrow(hc$merge)
  ind <- 1:(nr + 1)
  gr <- 1:(nr + 1)

  # for each merge (=level of hierarchical clustering)
  for (i in 1:nr)
  {
    # positive number = number of the line of the cluster merged
    # negative number = single variable

    indpos <- which(hc$merge[i, ] > 0)
    pos <- hc$merge[i, indpos]
    neg <- abs(hc$merge[i, ][hc$merge[i, ] < 0])

    # for each cluster, search the variables it contains
    while (length(pos) > 0) {
      pos <- hc$merge[pos, ]
      neg <- c(neg, abs(pos[pos < 0]))
      pos <- pos[pos > 0]
    }

    # complete output
    ind <- c(ind, sort(neg))
    gr <- c(gr, rep(i + nr + 1, length(neg)))
  } # fin for row merge

  return(list(varGroup = ind, indexGroup = gr))
}
