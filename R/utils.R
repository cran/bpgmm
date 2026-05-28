#' (internal)
#' @noRd
#' @noRd
get_z_mat_r <- function(ZOneDim, m, n) {
  Zmat <- matrix(NA, m, n)
  for (i in 1:m) {
    Zmat[i, ] <- as.numeric(ZOneDim == i)
  }
  Zmat
}

#' (internal)
#' @noRd
#' @noRd
calculate_ratio <- function(deno, nume) {
  ## deno nume both in log scale
  maxNume <- max(nume)
  transDeno <- deno - maxNume
  transNume <- nume - maxNume
  res <- exp(transDeno) / (sum(exp(transNume)))
  res
}


#' (internal)
#' @noRd
#' @noRd
constraint_list_to_models <- function(stringList) {
  vapply(stringList, constraint_to_model, character(1))
}



#' (internal)
#' @noRd
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

#' Summarize cluster allocations by posterior mode.
#'
#' @param Zlist List of sampled allocation vectors.
#' @param index Sample indices to include in the posterior mode.
#'
#' @return An integer vector of modal cluster allocations.
#' @noRd
#' @noRd
summarize_allocations <- function(Zlist, index = seq_along(Zlist)) {
  sampleSize <- length(Zlist[[1]])
  res <- integer(sampleSize)

  for (i in seq_len(sampleSize)) {
    temp <- integer(length(index))
    for (j in seq_along(index)) {
      temp[j] <- Zlist[[index[j]]][i]
    }
    res[i] <- get_mode(temp)
  }

  res
}

#' (internal)
#' @noRd
summarize_allocations_legacy <- summarize_allocations



#' (internal)
#' @noRd
#' @noRd
get_mode <- function(v) {
  v <- v[!is.na(v)]
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}


#' (internal)
#' @noRd
#' @noRd
calculate_covariance_list <- function(psyList, lambdaList) {
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
#' @noRd
update_latent_scores <- function(X, thetaYList, ZOneDim, clusInd, qVec) {
  update_latent_scores_cpp(X, thetaYList, ZOneDim, clusInd, qVec)
}


#' (internal)
#' @noRd
#' @noRd
clear_inactive_components <- function(thetaYList, clusInd, mMax) {
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
#' @noRd
combine_cluster_parameters <- function(oldList, newList, ind) {
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
#' @noRd
subset_theta_y <- function(thetaYList, Ind) {
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
#' @noRd
remove_theta_y_component <- function(thetaYList, Ind) {
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
#' @noRd
change_constraint_format <- function(strNum) {
  digits <- regmatches(strNum, gregexpr("[01]", strNum))[[1]]
  constraint_to_model(as.integer(digits))
}
