#' PriorThetaY list
#' @description generate prior value for parameter Theta and Y.
#' @importFrom gtools rdirichlet
#' @importFrom mvtnorm rmvnorm
#' @import stats
#' @param m the number of cluster
#' @param n sample size
#' @param p number of covariates
#' @param muBar parameter
#' @param hparam hyperparameters
#' @param qVec the vector of the number of factors in each clusters
#' @param ZOneDim ZOneDim
#' @param constraint constraint
#'
#' @return PriorThetaY list
#'
#' @examples
#'
#' sample_data = "https://raw.githubusercontent.com/lzyacht/bpgmm-examples/master/data/sampleData.csv"
#' X = utils::read.table(sample_data, header = TRUE, sep = ',')
#' X = as.matrix(X)
#' nsim = 1
#' burn = 20
#' n = ncol(X)
#' p = nrow(X)
#' m = 2
#' qVec = rep(3,m)
#' qnew = 3
#' delta = 2
#' ggamma = 2
#' dVec = c(1,1,1)
#' sVec = c(1,1,1)
#' constraint = c(1,1,1)
#' hparam = new("Hparam", alpha1=1, alpha2=2, bbeta=3, delta=4, ggamma=5)
#' muBar = apply(X, MARGIN = 1, FUN = mean)
#' ZOneDim = kmeans(x = t(X), centers = m)$cluster
#' generatePriorThetaY(m, n, p, muBar, hparam, qVec, ZOneDim, constraint)
#'
#'
#' @export
#'
generatePriorThetaY = function(m,
                               n,
                               p,
                               muBar,
                               hparam,
                               qVec,
                               ZOneDim,
                               constraint) {
  ggamma <- hparam@ggamma
  delta  <- hparam@delta
  bbeta  <- hparam@bbeta
  alpha1 <- hparam@alpha1
  alpha2 <- hparam@alpha2

  # prior tao
  if (m == 1) {
    tao <- rbeta(n = 1, shape1 = 1, shape2 = m)
  } else{
    tao <- gtools::rdirichlet(n = 1, alpha = rep(ggamma, m))
  }

  # prior psy
  psy = generatePriorPsi(p, m, delta, bbeta, constraint)

  # prior M
  M = list()
  for (i in 1:m) {
    M[[i]] = (mvtnorm::rmvnorm(n =  1, mean = muBar, sigma = 1 / alpha1 * psy[[i]]))
  }

  # prior lambda
  lambda = generatePriorLambda(p, m, alpha2, qVec, psy, constraint)


  # Zmat = getZmat(ZOneDim, m, n)
  Zmat = get_Z_mat(ZOneDim, m, n)

  # post Y
  Y = list()
  for (k in 1:m) {
    Y[[k]] = matrix(NA, qVec[k], n)
    for (i in 1:n) {
      Y[[k]][, i] = mvtnorm::rmvnorm(1, mean = rep(0, qVec[k]), sigma = diag(qVec[k]))
    }
  }

  new("ThetaYList",
      tao    = tao,
      psy    = psy,
      M      = M,
      lambda = lambda,
      Y      = Y)
}

#' evaluate Prior
#' @description evaluate prior value for parameter Theta and Y.
#' @importFrom gtools ddirichlet
#' @importFrom mvtnorm dmvnorm
#' @import stats
#'
#' @param m m
#' @param p p
#' @param muBar mu_bar
#' @param hparam hyper parameter class
#' @param thetaYList theta Y list
#' @param ZOneDim one dime of z
#' @param qVec q vector
#' @param constraint type of constraint
#'
#' @return evaluated prior value for parameter Theta and Y.
#'
#' @examples
#' \donttest{
#' url <- paste0("https://github.com/lzyacht/bpgmm-examples/",
#' "blob/master/data/updatePostThetaY_example.RData?raw=true")
#' download.file(url, destfile= "updatePostThetaY_example.RData", mode = "wb")
#' load("updatePostThetaY_example.RData")
#' evaluatePrior(m, p, muBar,hparam, thetaYList, ZOneDim, qVec,constraint)
#' }
#'
#' @export
evaluatePrior = function(m,
                         p,
                         muBar,
                         hparam,
                         thetaYList,
                         ZOneDim,
                         qVec,
                         constraint){

  ggamma = hparam@ggamma
  alpha1 = hparam@alpha1
  alpha2 = hparam@alpha2
  bbeta  = hparam@bbeta
  delta  = hparam@delta
  ## 2.2: Y
  Yval = 0
  for(k in 1:m){
    Yval = Yval + sum(dnorm(thetaYList@Y[[k]], log = T))
  }
  ## 2.3: Z
  adjustTao = thetaYList@tao/sum(thetaYList@tao)
  Zval = log(sum(adjustTao[ZOneDim]))
  ## 2.6: tao
  taoVal = log(gtools::ddirichlet(x = adjustTao, alpha = rep(ggamma, m)))

  ## 2.7: M
  Mval = 0
  for(k in 1:m){
    Mval = Mval + mvtnorm::dmvnorm(c(thetaYList@M[[k]]), mean = c(muBar), sigma = 1/alpha1 * thetaYList@psy[[k]], log = T)
  }

  ## 2.8: lambda
  lambdaVal = evaluatePriorLambda(p, m, alpha2, qVec,thetaYList@psy, thetaYList@lambda, constraint)

  ## 2.9: psy
  psyVal = evaluatePriorPsi(thetaYList@psy, p, m, delta, bbeta, constraint)

  totalVal = sum(Yval + Zval + taoVal + Mval + lambdaVal + psyVal)
  return(totalVal)

}
#' generatePriorPsi
#'
#' @description generate prior value for parameter Psi
#' @import stats
#' @param p the number of features
#' @param m the number of clusters
#' @param delta hyperparameters
#' @param bbeta hyperparameters
#' @param constraint the pgmm constraint, a vector of length three with binary entry. For example, c(1,1,1) means the fully constraint model
#' @return generated prior value for parameter Psi
#'
#'
generatePriorPsi = function(p,
                            m,
                            delta,
                            bbeta,
                            constraint){
  psy = list()

  if(constraint[2] == T & constraint[3] == T){

    for(i in 1:m){
      if(i == 1){
        psyValue = 1/rgamma(1, shape = delta, rate = bbeta)
        psy[[i]] = diag(rep(psyValue, p))
      }else{
        psy[[i]] = psy[[1]]
      }
    }

  }else if(constraint[2] == T & constraint[3] == F){
    for(i in 1:m){
      if(i == 1){
        psy[[i]] = diag(1/rgamma(p, shape = delta, rate = bbeta))
      }else{
        psy[[i]] =  psy[[1]]
      }
    }
  }else if(constraint[2] == F & constraint[3] == T){

    for(i in 1:m){
      psyValue = 1/rgamma(1, shape = delta, rate = bbeta)
      psy[[i]] = diag(rep(psyValue, p))
    }

  }else{
    for(i in 1:m){
      psy[[i]] = diag(1/rgamma(p, shape = delta, rate = bbeta))
    }
  }
  return(psy)
}

#' evaluatePriorPsi
#'
#' @description evaluate prior value for parameter Psi
#' @import stats
#' @param psy parameter
#' @param p the number of features
#' @param m the number of clusters
#' @param delta parameter
#' @param bbeta parameter
#' @param constraint parameter
#'
#' @return evaluated prior value for parameter Psi
#'
#'
evaluatePriorPsi = function(psy,
                            p,
                            m,
                            delta,
                            bbeta,
                            constraint){

  psyeval = 0
  if(constraint[2] == T & constraint[3] == T){

    for(i in 1:m){
      if(i == 1){
        psyValue = 1/psy[[i]][1,1]
        psyeval = psyeval + dgamma(psyValue, shape = delta, rate = bbeta, log = T)
      }
    }

  }else if(constraint[2] == T & constraint[3] == F){
    for(i in 1:m){
      if(i == 1){
        psyValue = 1/diag(psy[[i]])
        psyeval = psyeval + sum(dgamma(psyValue, shape = delta, rate = bbeta, log = T))
      }
    }
  }else if(constraint[2] == F & constraint[3] == T){

    for(i in 1:m){
      psyValue = 1/psy[[i]][1,1]
      psyeval = psyeval + dgamma(psyValue, shape = delta, rate = bbeta, log = T)
    }

  }else{
    for(i in 1:m){
      psyValue = 1/diag(psy[[i]])
      psyeval = psyeval + sum(dgamma(psyValue, shape = delta, rate = bbeta, log = T))
    }
  }
  return(psyeval)
}

#' generatePriorLambda
#'
#' @description evaluate prior value for parameter Lambda
#' @importFrom mvtnorm rmvnorm
#' @param p the number of features
#' @param m the number of clusters
#' @param alpha2 hyper parameter
#' @param qVec parameter
#' @param psy parameter
#' @param constraint parameter
#'
#' @return evaluated prior value for parameter Lambda
#'
generatePriorLambda = function(p,
                               m,
                               alpha2,
                               qVec,
                               psy,
                               constraint){

  lambda = list()
  if(constraint[1] == T & constraint[2] == T){
    for(k in 1:m){
      if(k == 1){
        qk = qVec[k]
        lambda[[k]] = matrix(0, p, qk)
        for(j in 1:qk){
          lambda[[k]][,j] = mvtnorm::rmvnorm(1, rep(0,p), 1/alpha2 * psy[[k]])
        }
      }else{
        lambda[[k]] = lambda[[1]]
      }
    }
  }else if(constraint[1] == T & constraint[2] == F){

    psyAve = matrix(0, nrow = p, ncol = p)
    for(k in 1:m){
      psyAve = psyAve +  solve(psy[[k]])
    }
    psyAve = solve(1/m * psyAve)

    for(k in 1:m){
      if(k == 1){
        qk = qVec[k]
        lambda[[k]] = matrix(0, p, qk)
        for(j in 1:qk){
          lambda[[k]][,j] = mvtnorm::rmvnorm(1, rep(0,p), 1/alpha2 * psyAve)
        }
      }else{
        lambda[[k]] = lambda[[1]]
      }
    }
  }else{
    for(k in 1:m){
      qk = qVec[k]
      lambda[[k]] = matrix(0, p, qk)
      for(j in 1:qk){
        lambda[[k]][,j] = mvtnorm::rmvnorm(1, rep(0,p), 1/alpha2 * psy[[k]])
      }
    }
  }
  return(lambda)
}

#' evaluatePriorLambda
#'
#' @description evaluate prior value for parameter Lambda
#' @importFrom mvtnorm dmvnorm
#' @param p the number of features
#' @param m the number of clusters
#' @param alpha2 hyper parameter
#' @param qVec the vector of the number of factors in each clusters
#' @param psy parameter
#' @param lambda parameter
#' @param constraint the pgmm constraint, a vector of length three with binary entry. For example, c(1,1,1) means the fully constraint model
#'
#' @return evaluated prior value for parameter Lambda
#'
evaluatePriorLambda = function(p,
                               m,
                               alpha2,
                               qVec,
                               psy,
                               lambda,
                               constraint){

  evallambda = 0
  if(constraint[1] == T & constraint[2] == T){
    for(k in 1:m){
      if(k == 1){
        qk = qVec[k]
        for(j in 1:qk){
          evallambda = evallambda + mvtnorm::dmvnorm(lambda[[k]][,j], rep(0,p), 1/alpha2 * psy[[k]], log = T)
        }
      }
    }
  }else if(constraint[1] == T & constraint[2] == F){

    psyAve = matrix(0, nrow = p, ncol = p)
    for(k in 1:m){
      psyAve = psyAve +  solve(psy[[k]])
    }
    psyAve = solve(1/m * psyAve)

    for(k in 1:m){
      if(k == 1){
        qk = qVec[k]
        for(j in 1:qk){
          evallambda = evallambda + mvtnorm::dmvnorm(lambda[[k]][,j], rep(0,p), 1/alpha2 * psyAve, log = T)
        }
      }
    }
  }else{
    for(k in 1:m){
      qk = qVec[k]
      for(j in 1:qk){
        evallambda = evallambda + mvtnorm::dmvnorm(lambda[[k]][,j], rep(0,p), 1/alpha2 * psy[[k]], log = T)
      }
    }
  }
  return(evallambda)
}
