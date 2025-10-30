#' (internal)
#' @noRd
VstepRJMCMCupdate <- function(X,
                              muBar,
                              p,
                              thetaYList,
                              ZOneDim,
                              hparam,
                              hparamInit,
                              qVec,
                              qnew,
                              ggamma,
                              dVec,
                              sVec,
                              constraint,
                              clusInd) {



  ## rjmcmc for V
  # print("vStep")
  n <- ncol(X)
  m <- sum(clusInd)
  # steps = c("stay", "lambda", "psi1", "psi2")
  steps <- c("lambda", "psi1", "psi2")
  currentStep <- sample(size = 1, x = steps)
  # currentStep = "psi2"
  # currentStep = "stay"


  ## trans to non empty obj
  NEobj <- toNEthetaYlist(thetaYList, ZOneDim, qVec, clusInd)
  thetaYList <- NEobj$thetaYList
  ZOneDim <- NEobj$ZOneDim
  qVec <- NEobj$qVec
  # if(currentStep == "stay"){
  #   print("stay")
  #   ZOneDim = update_PostZ(X, m, n, thetaYList)
  #   # ZOneDim = updatePostZ(m, n, thetaYList)
  #   thetaYList = updatePostThetaY(m,  n,hparam, thetaYList, ZOneDim,  qVec,constraint, X, ggamma)
  #   hparam = update_Hyperparameter(m, p, qnew, hparam,thetaYList, dVec, sVec)
  #
  # }
  if (currentStep == "lambda") {

    # print("lambda step")
    proposeConstraint <- constraint
    proposeConstraint[1] <- (constraint[1] + 1) %% 2

    CxyList <- Calculate_Cxy(m, n, hparam, thetaYList, ZOneDim, qVec, X)

    newLambda <- CalculateProposalLambda(hparam, thetaYList, CxyList, proposeConstraint, m, p, qVec)

    newthetaYList <- thetaYList
    newthetaYList@lambda <- newLambda

    newCxyList <- Calculate_Cxy(m, n, hparam, newthetaYList, ZOneDim, qVec, X)

    oldDensity <- likelihood(thetaYList, ZOneDim, qVec, muBar, X) + evaluatePrior(m, p, muBar, hparam, thetaYList, ZOneDim, qVec, constraint, sort(clusInd, decreasing = T))
    oldLambdaEval <- EvaluateProposalLambda(hparam, newthetaYList, newCxyList, constraint, thetaYList@lambda, m, qVec, p)
    newLambdaEval <- EvaluateProposalLambda(hparam, thetaYList, CxyList, proposeConstraint, newthetaYList@lambda, m, qVec, p)




    ## Gibbs
    newhparam <- hparam
    for (i in 1:10) {
      newZOneDim <- update_PostZ(X, m, n, newthetaYList)
      # newZOneDim = updatePostZ(m, n, newthetaYList)

      newthetaYList <- updatePostThetaY(m = m, n = n, p, newhparam, newthetaYList, ZOneDim = newZOneDim, qVec = qVec, constraint = proposeConstraint, X, ggamma)
      newhparam <- update_Hyperparameter(m, p, qnew, newhparam, newthetaYList, dVec, sVec)
    }

    newDensity <- likelihood(newthetaYList, newZOneDim, qVec, muBar, X) + evaluatePrior(m, p, muBar, newhparam, newthetaYList, newZOneDim, qVec, proposeConstraint, sort(clusInd, decreasing = T))

    numer <- newDensity + oldLambdaEval
    denom <- oldDensity + newLambdaEval
    # diffVec[h] = numer - denom
    # cat("accepting prob = exp ", numer - denom, " \n")

    probAlpha <- calculateRatio(numer, denom)
    acceptP <- min(1, probAlpha)
    res <- rbinom(1, size = 1, prob = acceptP)

    if (res == 1) {
      # print("lambda success=====>")
      # cat("lambda prob = ", probAlpha, "====>\n")

      constraint <- proposeConstraint
      thetaYList <- newthetaYList
      ZOneDim <- newZOneDim
      hparam <- newhparam
    } else {
      # print("lambda fail=====>")
    }
  }
  else if (currentStep == "psi1" | currentStep == "psi2") {

    # print(currentStep)
    proposeConstraint <- constraint
    if (currentStep == "psi1") {
      proposeConstraint[2] <- (constraint[2] + 1) %% 2
    } else {
      proposeConstraint[3] <- (constraint[3] + 1) %% 2
    }

    CxyList <- Calculate_Cxy(m, n, hparam, thetaYList, ZOneDim, qVec, X)
    newPsy <- CalculateProposalPsy(hparam, thetaYList, CxyList, proposeConstraint, m, p, qVec)

    newthetaYList <- thetaYList
    newthetaYList@psy <- newPsy
    newCxyList <- Calculate_Cxy(m, n, hparam, newthetaYList, ZOneDim, qVec, X)


    # oldDensity = likelihood(thetaYList, ZOneDim,qVec,muBar, X) + evaluatePrior(m, p, muBar, hparam, thetaYList, ZOneDim, qVec, constraint,clusInd)
    # newDensity = likelihood(newthetaYList, ZOneDim,qVec,muBar, X) + evaluatePrior(m, p, muBar,hparam, newthetaYList, ZOneDim, qVec, proposeConstraint,clusInd)

    oldDensity <- likelihood(thetaYList, ZOneDim, qVec, muBar, X) + evaluatePrior(m, p, muBar, hparam, thetaYList, ZOneDim, qVec, constraint, sort(clusInd, decreasing = T))
    newDensity <- likelihood(newthetaYList, ZOneDim, qVec, muBar, X) + evaluatePrior(m, p, muBar, hparam, newthetaYList, ZOneDim, qVec, proposeConstraint, sort(clusInd, decreasing = T))


    oldPsyEval <- EvaluateProposalPsy(hparam, newthetaYList, newCxyList, constraint, thetaYList@psy, m, p, qVec)
    newPsyEval <- EvaluateProposalPsy(hparam, thetaYList, CxyList, proposeConstraint, newthetaYList@psy, m, p, qVec)

    numer <- newDensity + oldPsyEval
    denom <- oldDensity + newPsyEval
    # print(numer - denom)
    # diffVec[h] = numer - denom

    probAlpha <- calculateRatio(numer, denom)
    acceptP <- min(1, probAlpha)
    res <- rbinom(1, size = 1, prob = acceptP)

    if (res == 1) {
      # print("psi success=====>")
      # cat("psi prob = ", probAlpha, "====>\n")
      constraint <- proposeConstraint
      thetaYList <- newthetaYList
      # psy = newPsy
    } else {

      # print("psi fail=====>")
    }
  }

  ## trans to empty obj
  Eobj <- toEthetaYlist(thetaYList, ZOneDim, qnew, clusInd)
  thetaYList <- Eobj$thetaYList
  ZOneDim <- Eobj$ZOneDim
  qVec <- Eobj$qVec



  return(list(thetaYList = thetaYList, ZOneDim = ZOneDim, hparam = hparam, constraint = constraint))
}
