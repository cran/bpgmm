#' stayMCMCupdate
#'
#' @param X X
#' @param thetaYList thetaYList
#' @param ZOneDim ZOneDim
#' @param hparam hparam
#' @param qVec qVec
#' @param qnew qnew
#' @param dVec dVec
#' @param sVec sVec
#' @param constraint constraint
#' @param clusInd clusInd
#'
stayMCMCupdate <- function(X,
                           thetaYList,
                           ZOneDim,
                           hparam,
                           qVec,
                           qnew,
                           dVec,
                           sVec,
                           constraint,
                           clusInd) {
  ggamma <- hparam@ggamma
  m <- length(which(clusInd == 1))
  p <- nrow(X)
  n <- ncol(X)
  ## transfor to non empty obj
  sampleVec <- which(clusInd == 1)
  NEthetaYList <- getIndThetaY(thetaYList, sampleVec)

  NEZOneDim <- ZOneDim
  for (i in 1:length(sampleVec)) {
    NEZOneDim[NEZOneDim == sampleVec[i]] <- i
  }
  ## update

  NEthetaYList <- updatePostThetaY(m, n, p, hparam, NEthetaYList, NEZOneDim, qVec[qVec != 0], constraint, X, ggamma)
  NEZOneDim <- updatePostZ(X, m, n, NEthetaYList)
  hparam <- update_Hyperparameter(m, p, qnew, hparam, NEthetaYList, dVec, sVec)


  ZOneDim <- sampleVec[NEZOneDim]

  thetaYList <- getThetaYWithEmpty(NEthetaYList, clusInd)
  return(list(thetaYList = thetaYList, ZOneDim = ZOneDim, hparam = hparam))
}



#' (internal)
#' @noRd
getThetaYWithEmpty <- function(NEthetaYList, clusInd) {
  thetaYList <- new("ThetaYList")
  j <- 1
  for (i in 1:length(clusInd)) {
    if (clusInd[i] == 1) {
      thetaYList@tao[i] <- NEthetaYList@tao[j]
      thetaYList@psy[[i]] <- NEthetaYList@psy[[j]]
      thetaYList@M[[i]] <- NEthetaYList@M[[j]]
      thetaYList@lambda[[i]] <- NEthetaYList@lambda[[j]]
      thetaYList@Y[[i]] <- NEthetaYList@Y[[j]]
      j <- j + 1
    } else {
      thetaYList@tao[i] <- 0
      thetaYList@psy[[i]] <- NA
      thetaYList@M[[i]] <- NA
      thetaYList@lambda[[i]] <- NA
      thetaYList@Y[[i]] <- NA
    }
  }
  thetaYList
}

#' (internal)
#' @noRd
toNEthetaYlist <- function(thetaYList, ZOneDim, qVec, clusInd) {
  ## transfor to non empty obj
  sampleVec <- which(clusInd == 1)
  NEthetaYList <- getIndThetaY(thetaYList, sampleVec)

  NEZOneDim <- ZOneDim
  for (i in 1:length(sampleVec)) {
    NEZOneDim[NEZOneDim == sampleVec[i]] <- i
  }
  qnew <- max(qVec)
  qVec <- qVec[qVec == qnew]
  return(list(thetaYList = NEthetaYList, ZOneDim = NEZOneDim, qVec = qVec))
}

#' Title
#'
#' @param NEthetaYList NEthetaYList
#' @param NEZOneDim NEZOneDim
#' @param qnew qnew
#' @param clusInd clusInd
#'
toEthetaYlist <- function(NEthetaYList, NEZOneDim, qnew, clusInd) {
  ## transfor back
  sampleVec <- which(clusInd == 1)
  ZOneDim <- NEZOneDim
  for (i in length(sampleVec):1) {
    ZOneDim[ZOneDim == i] <- sampleVec[i]
  }
  qVec <- clusInd
  qVec[clusInd == 1] <- qnew
  thetaYList <- getThetaYWithEmpty(NEthetaYList, clusInd)
  return(list(thetaYList = thetaYList, ZOneDim = ZOneDim, qVec = qVec))
}
