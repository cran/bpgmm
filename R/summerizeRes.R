#' summerizePgmmRJMCMC
#'
#' @param pgmmResList result list from pgmmRJMCMC
#' @param trueCluster true cluster allocation
#' @importFrom mclust adjustedRandIndex
#'
summerizePgmmRJMCMC <- function(pgmmResList, trueCluster = NULL) {
  Zalloc <- sumerizeZ(res$ZmatList)

  nCluster <- table(sapply(res$ZmatList, function(x) {
    length(unique(x))
  }))

  nConstraint <- res$constraintList
  nConstraint <- listToStrVec(nConstraint)
  nConstraint <- table(nConstraint, dnn = "")



  sumRes <- list(Zalloc = Zalloc, nCluster = nCluster, nConstraint = nConstraint)

  if (!is.null(trueCluster)) {
    ari <- adjustedRandIndex(trueCluster, Zalloc)
    sumRes$ari <- ari
  }

  sumRes
}
