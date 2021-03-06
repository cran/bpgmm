#' Update posterior theta Y list
#'
#' @param m the number of clusters.
#' @param p the number of variables
#' @param n the number of observations.
#' @param hparam hyper parameters
#' @param thetaYList theta Y list
#' @param ZOneDim ZOneDim
#' @param qVec qVec
#' @param constraint constraint
#' @param X X
#' @param ggamma ggamma
#'
#' @examples
#' set.seed(100)
#' n <- 10
#' p <- 2
#' q <- 1
#' K <- 2
#' m <- 1
#' muBar <- c(0, 0)
#' qVec <- c(1, 1)
#' constraint <- c(0, 0, 0)
#' X <- t(
#'   fabMix::simData(
#'     sameLambda = TRUE,
#'     sameSigma = TRUE,
#'     K.true = K,
#'     n = n,
#'     q = q,
#'     p = p,
#'     sINV_values = 1 / ((1:p))
#'   )$data
#' )
#' hparam <- new(
#'   "Hparam",
#'   alpha1 = 0.567755037123148,
#'   alpha2 = 1.1870201935945,
#'   delta = 2,
#'   ggamma = 2,
#'   bbeta = 3.39466184520673
#' )
#' ZOneDim <- sample(seq_len(m), n, replace = TRUE)
#' thetaYList <-
#'   new(
#'     "ThetaYList",
#'     tao = 0.366618687752634,
#'     psy = list(structure(
#'       c(
#'         4.18375613018654,
#'         0, 0, 5.46215996830771
#'       ),
#'       .Dim = c(2L, 2L)
#'     )),
#'     M = list(structure(
#'       c(
#'         3.27412045866392,
#'         -2.40544145363349
#'       ),
#'       .Dim = 1:2
#'     )),
#'     lambda = list(structure(
#'       c(
#'         2.51015961514781,
#'         -0.0741189919182549
#'       ),
#'       .Dim = 2:1
#'     )),
#'     Y = list(structure(
#'       c(
#'         -0.244239011725104,
#'         -0.26876172736886,
#'         0.193431511203083,
#'         0.41624466812811,
#'         -0.54581548068437,
#'         -0.0479517628308146,
#'         -0.633383997203325,
#'         0.856855296613208,
#'         0.792850576988512,
#'         0.268208848994559
#'       ),
#'       .Dim = c(1L, 10L)
#'     ))
#'   )
#' constraint <- c(0, 0, 0)
#' #'
#' \donttest{
#' updatePostThetaY(m, n, p, hparam, thetaYList, ZOneDim, qVec, constraint, X, ggamma)
#' }
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
