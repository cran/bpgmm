#' (internal)
#' @noRd
updatePostThetaY <- function(m, n, p, hparam, thetaYList, ZOneDim, qVec, constraint, X, ggamma) {
  alpha1 <- hparam@alpha1
  alpha2 <- hparam@alpha2
  bbeta <- hparam@bbeta
  lambda <- thetaYList@lambda
  Y <- thetaYList@Y
  M <- thetaYList@M
  psy <- thetaYList@psy

  ## post for Theta = {tao, M, Lambda, psy}
  CxyList <- Calculate_Cxy(m, n, hparam, thetaYList, ZOneDim, qVec, X)
  Cxxk <- CxyList$Cxxk
  Cxyk <- CxyList$Cxyk
  Cyyk <- CxyList$Cyyk
  Cytytk <- CxyList$Cytytk
  Cxtytk <- CxyList$Cxtytk
  CxL1k <- CxyList$CxL1k
  Cxmyk <- CxyList$Cxmyk
  sumCxmyk <- CxyList$sumCxmyk
  sumCyyk <- CxyList$sumCyyk
  A <- CxyList$A
  nVec <- CxyList$nVec

  Zmat <- getZmat(ZOneDim, m, n)


  # post tao
  tao <- rdirichlet(1, nVec + ggamma)

  #  post mu
  M <- list()
  for (k in 1:m) {
    M[[k]] <- rmvnorm(1, mean = CxL1k[[k]] * (sum(Zmat[k, ]) + alpha1)^(-1), sigma = (sum(Zmat[k, ]) + alpha1)^(-1) * psy[[k]])
  }

  ## lambda; psy
  lambdaPsy <- Calculate_PostLambdaPsy(m, p, hparam, CxyList, thetaYList, qVec, constraint)
  # lambdaPsy = CalculatePostLambdaPsy(alpha1, alpha2, bbeta, CxyList, M, psy, constraint)
  lambda <- lambdaPsy$lambda
  psy <- lambdaPsy$psy

  ## post Y
  D <- list()
  for (i in 1:m) {
    D[[i]] <- t(lambda[[i]]) %*% solve(psy[[i]] + lambda[[i]] %*% t(lambda[[i]]), tol = 1e-100)
  }

  Sigma <- list()
  for (i in 1:m) {
    Sigma[[i]] <- diag(qVec[i]) - D[[i]] %*% lambda[[i]]
  }

  Y <- list()
  YDvalList <- c()
  for (k in 1:m) {
    Y[[k]] <- matrix(NA, qVec[k], n)
    for (i in 1:n) {
      if (Zmat[k, i] == 0) {
        Y[[k]][, i] <- rmvnorm(1, mean = rep(0, qVec[k]), sigma = diag(qVec[k]))
      } else if (Zmat[k, i] == 1) {
        Y[[k]][, i] <- rmvnorm(1, mean = D[[k]] %*% t(X[, i] - M[[k]]), sigma = Sigma[[k]])
      }
    }
  }

  new("ThetaYList", tao = tao, psy = psy, M = M, lambda = lambda, Y = Y)
}
