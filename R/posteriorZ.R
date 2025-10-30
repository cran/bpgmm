#' (internal)
#' @noRd
updatePostZ <- function(X, m, n, thetaYList) {
  tao <- thetaYList@tao
  psy <- thetaYList@psy
  M <- thetaYList@M
  lambda <- thetaYList@lambda


  pMat <- matrix(NA, m, n)
  ## evaluate density
  dMat <- matrix(NA, m, n)

  for (k in 1:m) {
    for (i in 1:n) {
      dMat[k, i] <- mvtnorm::dmvnorm(X[, i],
        mean = M[[k]], sigma = psy[[k]] + lambda[[k]] %*% t(lambda[[k]]),
        log = T
      )
    }
  }

  for (k in 1:m) {
    dMat[k, ] <- dMat[k, ] + log(tao[k])
  }

  for (i in 1:n) {
    for (k in 1:m) {
      pMat[k, i] <- calculateRatio(dMat[k, i], dMat[, i])
    }
  }

  ZOneDim <- c()
  for (i in 1:n) {
    tempProb <- as.numeric(pMat[, i])
    ZOneDim[i] <- sample(x = 1:m, size = 1, prob = tempProb)
  }
  ZOneDim
}
