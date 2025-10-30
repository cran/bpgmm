#' (internal)
#' @noRd
getZmat <- function(ZOneDim, m, n) {
  Zmat <- matrix(NA, m, n)
  for (i in 1:m) {
    Zmat[i, ] <- as.numeric(ZOneDim == i)
  }
  Zmat
}

#' (internal)
#' @noRd
calculateRatio <- function(deno, nume) {
  ## deno nume both in log scale
  maxNume <- max(nume)
  transDeno <- deno - maxNume
  transNume <- nume - maxNume
  res <- exp(transDeno) / (sum(exp(transNume)))
  res
}


#' (internal)
#' @noRd
listToStrVec <- function(stringList) {
  for (i in 1:length(stringList)) {
    stringList[[i]] <-
      paste0("(", paste0(stringList[[i]], collapse = ""), ")")
  }
  unlist(stringList)
}



#' (internal)
#' @noRd
likelihood <- function(thetaYList, ZOneDim, qqVec, muBar, X) {
  m <- length(qqVec)
  n <- ncol(X)
  xvec <- c()
  ## 2.1: X
  xVal <- 0
  for (i in 1:n) {
    ## for subject i
    k <- ZOneDim[i]

    # if(is.vector(X[, ZOneDim == k])){
    #   meanx = X[, ZOneDim == k]
    # }else{
    #   meanx = apply(X[, ZOneDim == k], 1, mean)
    # }

    meanx <- thetaYList@lambda[[k]] %*% thetaYList@Y[[k]][, i] + c(thetaYList@M[[k]])
    varx <- thetaYList@psy[[k]]
    xvec[i] <- (mvtnorm::dmvnorm(x = X[, i], mean = meanx, sigma = varx, log = T))
    xVal <- xVal + mvtnorm::dmvnorm(x = X[, i], mean = meanx, sigma = varx, log = T)
  }

  xVal
  # xvec
}

#' (internal)
#' @noRd
sumerizeZ <- function(Zlist, index = 1:length(Zlist)) {
  sampleSize <- length(Zlist[[1]])
  res <- c()
  for (i in 1:sampleSize) {
    temp <- c()
    for (j in index) {
      temp[j] <- Zlist[[j]][i]
    }
    res[i] <- getmode(temp)
  }
  return(res)
}



#' (internal)
#' @noRd
getmode <- function(v) {
  v <- v[!is.na(v)]
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}


#' (internal)
#' @noRd
calculateVarList <- function(psyList, lambdaList) {
  m <- length(psyList)
  varList <- list()
  for (i in 1:m) {
    tempCov <- psyList[[i]] + lambdaList[[i]] %*% t(lambdaList[[i]])
    varList[[i]] <- tempCov
  }
  return(varList)
}


#' (internal)
#' @noRd
clearCurrentThetaYlist <- function(thetaYList, clusInd, mMax) {
  resThetaYList <- thetaYList

  for (i in 1:mMax) {
    if (clusInd[i] == 1) {
      resThetaYList@tao[i] <- thetaYList@tao[i]
      resThetaYList@psy[[i]] <- thetaYList@psy[[i]]
      resThetaYList@M[[i]] <- thetaYList@M[[i]]
      resThetaYList@lambda[[i]] <- thetaYList@lambda[[i]]
      resThetaYList@Y[[i]] <- thetaYList@Y[[i]]
    } else {
      resThetaYList@tao[i] <- 0
      resThetaYList@psy[[i]] <- NA
      resThetaYList@M[[i]] <- NA
      resThetaYList@lambda[[i]] <- NA
      resThetaYList@Y[[i]] <- NA
    }
  }
  return(resThetaYList)
}


#' (internal)
#' @noRd
combineClusterPara <- function(oldList, newList, ind) {
  resList <- oldList
  resList@tao <- oldList@tao * (1 - newList@tao)
  resList@tao[ind] <- newList@tao

  resList@psy[ind] <- newList@psy
  resList@M[ind] <- newList@M
  resList@lambda[ind] <- newList@lambda
  resList@Y[ind] <- newList@Y
  return(resList)
}




#' (internal)
#' @noRd
getIndThetaY <- function(thetaYList, Ind) {
  new("ThetaYList",
    tao = thetaYList@tao[Ind],
    psy = thetaYList@psy[Ind],
    M = thetaYList@M[Ind],
    lambda = thetaYList@lambda[Ind],
    Y = thetaYList@Y[Ind]
  )
}

#' (internal)
#' @noRd
getRemovedIndThetaY <- function(thetaYList, Ind) {
  res <- thetaYList

  res@tao[Ind] <- NA
  res@tao <- res@tao / sum(res@tao, na.rm = T)
  res@psy[[Ind]] <- NA
  res@M[[Ind]] <- NA
  res@lambda[[Ind]] <- NA
  res@Y[[Ind]] <- NA

  res
}

#' (internal)
#' @noRd
changeConstraintFormat <- function(strNum) {
  res <- ""
  res <- gsub(pattern = "0", replacement = "U", x = strNum)
  res <- gsub(pattern = "1", replacement = "C", x = res)
  res <- gsub(pattern = "\\(", replacement = "", x = res)
  res <- gsub(pattern = "\\)", replacement = "", x = res)
  return(res)
}
