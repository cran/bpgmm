#' (internal)
#' @noRd
#' @noRd
update_post_theta_y <- function(m, n, p, hparam, thetaYList, ZOneDim, qVec, constraint, X, ggamma) {
  alpha1 <- hparam@alpha1
  alpha2 <- hparam@alpha2
  bbeta <- hparam@bbeta
  lambda <- thetaYList@lambda
  Y <- thetaYList@Y
  M <- thetaYList@M
  psy <- thetaYList@psy

  ## post for Theta = {tao, M, Lambda, psy}
  CxyList <- calculate_cxy(m, n, hparam, thetaYList, ZOneDim, qVec, X)
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

  Zmat <- get_z_mat_r(ZOneDim, m, n)


  # post tao
  tao <- rdirichlet(1, nVec + ggamma)

  #  post mu
  M <- list()
  for (k in 1:m) {
    M[[k]] <- rmvnorm(1, mean = CxL1k[[k]] * (sum(Zmat[k, ]) + alpha1)^(-1), sigma = (sum(Zmat[k, ]) + alpha1)^(-1) * psy[[k]])
  }

  ## lambda; psy
  lambdaPsy <- calculate_post_lambda_psi(m, p, hparam, CxyList, thetaYList, qVec, constraint)
  # lambdaPsy = CalculatePostLambdaPsy(alpha1, alpha2, bbeta, CxyList, M, psy, constraint)
  lambda <- lambdaPsy$lambda
  psy <- lambdaPsy$psy

  ## post Y
  theta_for_y <- new("ThetaYList", tao = tao, psy = psy, M = M, lambda = lambda, Y = Y)
  Y <- update_latent_scores_cpp(X, theta_for_y, ZOneDim, rep(1, m), qVec)

  new("ThetaYList", tao = tao, psy = psy, M = M, lambda = lambda, Y = Y)
}
