#' bpgmm Model-Based Clustering Using Baysian PGMM Carries out model-based clustering using parsimonious Gaussian mixture models. MCMC are used for parameter estimation. The RJMCMC is used for model selection.
#'
#' @import stats MASS mcmcse pgmm label.switching
#' @param niter the number of iterations
#' @param burn the number of burn in iterations
#' @param X the observation matrix with size p * m
#' @param n the number of observations
#' @param p the number of features
#' @param m the number of clusters
#' @param delta scaler hyperparameters
#' @param ggamma scaler hyperparameters
#' @param qVec the vector of the number of factors in each clusters
#' @param qnew the number of factor for a new cluster
#' @param constraint the pgmm constraint, a vector of length three with binary entry. For example, c(1,1,1) means the fully constraint model
#' @param dVec a vector of hyperparameters with length three, shape parameters for alpha1, alpha2 and bbeta respectively
#' @param sVec a vector of hyperparameters with length three, rate parameters for alpha1, alpha2 and bbeta respectively
#'
#' @return parsimonious Gaussian mixture models classification results list
#'
#' @examples
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
#' constraint = c(0,0,0)
#' \donttest{
#' parsimoniousGaussianMixtureModel(nsim,burn,X,n,p,delta,ggamma,m,qVec,qnew,constraint,dVec,sVec)
#' }
#'
#'
#' @export
parsimoniousGaussianMixtureModel = function(niter,
                                            burn,
                                            X,
                                            n,
                                            p,
                                            delta,
                                            ggamma,
                                            m,
                                            qVec,
                                            qnew,
                                            constraint,
                                            dVec,
                                            sVec){


  alpha1 = rgamma(1, dVec[1], sVec[1])
  alpha2 = rgamma(1, dVec[2], sVec[2])
  bbeta  = rgamma(1, dVec[3], sVec[3])

  hparam = new("Hparam", alpha1=alpha1, alpha2=alpha2, bbeta=bbeta, delta=delta, ggamma=ggamma)

  muBar = apply(X, MARGIN = 1, FUN = mean)

  ## priors
  ZOneDim = kmeans(x = t(X), centers = m)$cluster
  thetaYList = generatePriorThetaY(m, n, p, muBar, hparam, qVec, ZOneDim, constraint)

  # thetaYList = new("ThetaYList",tao = priorList$tao,
  #                  psy = priorList$psy, M = priorList$M,
  #                  lambda = priorList$lambda, Y = priorList$Y
  #                  )

  for(i in 1:burn){
    ZOneDim = update_PostZ(X, m, n, thetaYList)
    thetaYList = updatePostThetaY(m, n, p, hparam, thetaYList, ZOneDim,  qVec,constraint, X)
    hparam = update_Hyperparameter(m, p, qnew, hparam,thetaYList, dVec, sVec)
  }

  ##
  alpha1Vec = c()
  alpha2Vec = c()
  bbetaVec = c()
  taoList = list()
  psyList = list()
  MList = list()
  lambdaList = list()
  YList = list()
  ZmatList = list()
  constraintList = list()
  ##
  diffVec = c()
  lambdaCount = 0

  ## posteriors
  for(h in 1:niter){
    cat("h = ", h, "  =========>\n")
    steps = c("stay", "lambda", "psi1", "psi2")
    currentStep = sample(size = 1, x = steps)
    # currentStep = "psi1"
    if(currentStep == "stay"){
      print("stay")
      ZOneDim = update_PostZ(X, m, n, thetaYList)
      thetaYList = updatePostThetaY(m, n, p, hparam, thetaYList, ZOneDim,  qVec,constraint, X)
      hparam = update_Hyperparameter(m, p, qnew, hparam,thetaYList, dVec, sVec)

    }else if(currentStep == "lambda"){

      print("lambda")
      proposeConstraint = constraint
      proposeConstraint[1] = (constraint[1] + 1) %% 2

      CxyList = Calculate_Cxy(m, n, hparam, thetaYList,ZOneDim, qVec, X)

      newLambda = CalculateProposalLambda(hparam, thetaYList, CxyList, proposeConstraint, m, p, qVec)

      newthetaYList = thetaYList
      newthetaYList@lambda = newLambda

      newCxyList = Calculate_Cxy(m, n, hparam,newthetaYList, ZOneDim, qVec, X)

      oldDensity = likelihood(thetaYList, ZOneDim,qVec,muBar, X) + evaluatePrior(m, p, muBar,hparam, thetaYList, ZOneDim, qVec,constraint)
      oldLambdaEval = EvaluateProposalLambda(hparam, newthetaYList,newCxyList, constraint, thetaYList@lambda, m, qVec, p)
      newLambdaEval = EvaluateProposalLambda(hparam, thetaYList,CxyList, proposeConstraint, newthetaYList@lambda, m, qVec, p)




      ## Gibbs
      newhparam = hparam
      for(i in 1:10){
        newZOneDim = update_PostZ(X, m, n, newthetaYList)
        newthetaYList = updatePostThetaY(m = m, n = n, p, newhparam, newthetaYList, ZOneDim=ZOneDim, qVec = qVec, constraint = proposeConstraint, X)
        newhparam = update_Hyperparameter(m, p, qnew, newhparam,newthetaYList, dVec, sVec)
      }

      newDensity = likelihood(newthetaYList, newZOneDim,qVec,muBar, X) + evaluatePrior(m, p, muBar, newhparam,newthetaYList, newZOneDim, qVec, proposeConstraint)

      numer  = newDensity + oldLambdaEval
      denom = oldDensity + newLambdaEval
      diffVec[h] = numer - denom
      cat("accepting prob = exp ", numer - denom, " \n")

      probAlpha = calculateRatio(numer, denom)
      acceptP = min(1,probAlpha)
      res = rbinom(1, size = 1, prob = acceptP)

      if(res == 1){
        lambdaCount = lambdaCount + 1
        print("lambda success=====>")
        cat("lambda prob = ", probAlpha, "====>\n")

        constraint = proposeConstraint
        thetaYList = newthetaYList
        ZOneDim = newZOneDim
        hparam = newhparam
      }
    }else if(currentStep == "psi1" | currentStep == "psi2"){

      print(currentStep)
      proposeConstraint = constraint
      if(currentStep == "psi1"){
        proposeConstraint[2] = (constraint[2] + 1) %% 2
      }else{
        proposeConstraint[3] = (constraint[3] + 1) %% 2
      }

      #####
      CxyList = Calculate_Cxy(m, n, hparam, thetaYList , ZOneDim, qVec, X)
      newPsy = CalculateProposalPsy(hparam, thetaYList, CxyList, proposeConstraint, m, p, qVec)

      newthetaYList = thetaYList
      newthetaYList@psy = newPsy
      newCxyList = Calculate_Cxy(m, n, hparam, newthetaYList, ZOneDim, qVec, X)


      oldDensity = likelihood(thetaYList, ZOneDim,qVec,muBar, X) + evaluatePrior(m, p, muBar, hparam, thetaYList, ZOneDim, qVec, constraint)
      newDensity = likelihood(newthetaYList, ZOneDim,qVec,muBar, X) + evaluatePrior(m, p, muBar,hparam, newthetaYList, ZOneDim, qVec, proposeConstraint)

      oldPsyEval = EvaluateProposalPsy(hparam,newthetaYList, newCxyList, constraint, thetaYList@psy, m, p, qVec)
      newPsyEval = EvaluateProposalPsy(hparam,thetaYList, CxyList, proposeConstraint, newthetaYList@psy,m, p, qVec)

      numer  = newDensity + oldPsyEval
      denom = oldDensity + newPsyEval
      print(numer - denom)
      diffVec[h] = numer - denom

      probAlpha = calculateRatio(numer, denom)
      acceptP = min(1,probAlpha)
      res = rbinom(1, size = 1, prob = acceptP)

    if(res == 1){
      print("psi success=====>")
      cat("psi prob = ", probAlpha, "====>\n")
      constraint = proposeConstraint
      psy = newPsy
    }
  }

    print(constraint)
    # print(ZOneDim)
    ##save
    alpha1Vec[h] = hparam@alpha1
    alpha2Vec[h] = hparam@alpha2
    bbetaVec[h] = hparam@bbeta
    taoList[[h]] = thetaYList@tao
    psyList[[h]] = thetaYList@psy
    MList[[h]] = thetaYList@M
    lambdaList[[h]] = thetaYList@lambda
    YList[[h]] = thetaYList@Y
    ZmatList[[h]] = ZOneDim
    constraintList[[h]] = constraint
  }

   switchedRed = lableSwichingByCluster(X,taoList,YList, MList, psyList,lambdaList, m, ZmatList)
   listTable = table(listToStrVec(constraintList))
   Zfinal = sumerizeZ(switchedRed$ZswitchList)
   return(list(listTable = listTable, Zfinal = Zfinal))

  # return(list(taoList=taoList, psyList=psyList,MList=MList,lambdaList=lambdaList,
  #          YList=YList,ZmatList=ZmatList, constraintList = constraintList,lambdaCount = lambdaCount,
  #            alpha1Vec = alpha1Vec, alpha2Vec = alpha2Vec, bbetaVec = bbetaVec, diffVec = diffVec))
}
