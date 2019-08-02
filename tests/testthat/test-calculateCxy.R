context("Utils - bpgmm calculateCxy")

test_that("calculateCxy is working correctly", {
  # R version of calculateCxy
  r_CalculateCxy <-  function(m, n, hparam, thetaYList, ZOneDim, qVec, X){

    alpha1 = hparam@alpha1
    alpha2 = hparam@alpha2
    bbeta  = hparam@bbeta

    Y      = thetaYList@Y
    lambda = thetaYList@lambda
    M      = thetaYList@M
    psy    = thetaYList@psy
    Zmat   = getZmat(ZOneDim,m,n)

    A      = list()
    nVec   = c()

    for(k in 1:m){
      nVec[k] = sum(Zmat[k,])
      A[[k]]  = diag(c(alpha1, alpha2 * rep(1, qVec[k])))
    }

    Cxxk = list()
    Cxyk = list()
    Cyyk = list()
    Cytytk = list()
    Cxtytk = list()
    CxL1k = list()
    Cxmyk = list()

    for(k in 1:m){
      Cxxk[[k]] = 0
      Cxyk[[k]] = 0
      Cyyk[[k]] = 0
      Cytytk[[k]] = 0
      CxL1k[[k]] = 0
      Cxmyk[[k]] = 0
      Cxtytk[[k]] = 0
      for(i in 1:n){
        Cxxk[[k]] = Cxxk[[k]] + Zmat[k,i] * X[,i] %*% t(X[,i])
        Cxyk[[k]] = Cxyk[[k]] + Zmat[k,i] * X[,i] %*% t(c(Y[[k]][,i]))
        Cyyk[[k]] = Cyyk[[k]] +  Zmat[k,i] * c(Y[[k]][,i]) %*% t(c(Y[[k]][,i]))
        Cxtytk[[k]] = Cxtytk[[k]] + Zmat[k,i] * X[,i] %*% t(c(1,Y[[k]][,i]))
        Cytytk[[k]] = Cytytk[[k]] + Zmat[k,i] * c(1,Y[[k]][,i]) %*% t(c(1,Y[[k]][,i]))

        Cxmyk[[k]] = Cxmyk[[k]] + Zmat[k,i] * c(X[,i] - M[[k]]) %*% t(c(Y[[k]][,i]))
        CxL1k[[k]] = CxL1k[[k]] + Zmat[k,i] * (X[,i] - lambda[[k]]%*%c(Y[[k]][,i]))
      }
    }

    sumCxmyk = 0
    sumCyyk = 0
    for(k in 1:m){
      sumCxmyk = sumCxmyk + Cxmyk[[k]]
      sumCyyk = sumCyyk + Cyyk[[k]] + alpha2*diag(qVec[k])
    }
    return(list(A = A,
                nVec = t(t(nVec)),
                Cxxk = Cxxk,
                Cxyk = Cxyk,
                Cyyk = Cyyk,
                Cytytk = Cytytk,
                Cxtytk = Cxtytk,
                CxL1k = CxL1k,
                Cxmyk = Cxmyk,
                sumCxmyk = sumCxmyk,
                sumCyyk = sumCyyk))
  }

  m = 3
  n = 7
  p = 10
  muBar =  rnorm(p)
  qVec = c(2,2,2)
  ZOneDim = c(1,2,3,1,2,2,3)
  constraint = c(1,1,1)
  hparam <- new("Hparam",alpha1 = 3, alpha2 = 2, delta  = 3, ggamma = 4, bbeta  = 5)
  thetaYList = generatePriorThetaY(m, n, p, muBar, hparam, qVec, ZOneDim, constraint)
  X = matrix(rnorm(p * n, 0, 1), p, n)
  res   = Calculate_Cxy(m, n, hparam, thetaYList, ZOneDim, qVec, X)
  r_res = r_CalculateCxy(m, n, hparam, thetaYList, ZOneDim, qVec, X)

  expect_equal(r_res$A, res$A)
  expect_equal(r_res$nVec, res$nVec)
  expect_equal(r_res$Cxxk, res$Cxxk)
  expect_equal(r_res$Cxyk, res$Cxyk)
  expect_equal(r_res$Cyyk, res$Cyyk)
  expect_equal(r_res$Cytytk, res$Cytytk)
  expect_equal(r_res$Cxtytk, res$Cxtytk)
  expect_equal(r_res$CxL1k, res$CxL1k)
  expect_equal(r_res$Cxmyk, res$Cxmyk)
  expect_equal(r_res$sumCxmyk, res$sumCxmyk)
  expect_equal(r_res$sumCyyk, res$sumCyyk)
})

