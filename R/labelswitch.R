library("label.switching")
lableSwichingByCluster = function(X,
                                  taoList,
                                  Ylist,
                                  Mlist,
                                  psyList,
                                  lambdaList,
                                  m,
                                  Zlist) {
  nsim = length(Mlist)
  probList = list()
  for (i in 1:nsim) {
    # m = mVec[[i]]
    distValue = matrix(NA, ncol(X), m)
    muList = Mlist[[i]]
    varList = calculateVarList(psyList[[i]], lambdaList[[i]])
    for (j in 1:m) {
      tempMean =  muList[[j]]
      tempVar = varList[[j]]
      funList = function(x) {
        dmvnorm(x = x,
                mean = tempMean,
                sigma =  tempVar)
      }
      distValue[, j] = apply(X = t(X),
                             MARGIN = 1,
                             FUN = funList)
    }
    probList[[i]] = distValue
  }

  ZswitchList = Zlist
  taoSwitchList = taoList
  MswitchList = Mlist
  YswitchList = Ylist
  lambdaSwitchList = lambdaList
  psySwitchList = psyList

  ind = 1:length(taoList)
  subList = probList[ind]
  pAry = array(dim = c(length(subList), ncol(X), m))
  for (j in 1:length(subList)) {
    pAry[j, , ] = subList[[j]]
  }

  if (m != 1) {
    perm = stephens(pAry)$permutations
    for (i in 1:length(ind)) {
      ZswitchList[[ind[i]]] = switchZFun(Zlist[[ind[i]]], perm[i, ])
      taoSwitchList[[ind[i]]] = switchTaoFun(taoList[[ind[i]]], perm[i, ])
      MswitchList[[ind[i]]] = swithListFun(Mlist[[ind[i]]], perm[i, ])
      YswitchList[[ind[i]]] = swithListFun(Ylist[[ind[i]]], perm[i, ])
      lambdaSwitchList[[ind[i]]] = swithListFun(lambdaList[[ind[i]]], perm[i, ])
      psySwitchList[[ind[i]]] = swithListFun(psyList[[ind[i]]], perm[i, ])
    }
  }

  return(
    list(
      ZswitchList = ZswitchList,
      taoSwitchList = taoSwitchList,
      MswitchList = MswitchList,
      YswitchList = YswitchList,
      lambdaSwitchList = lambdaSwitchList,
      psySwitchList = psySwitchList
    )
  )
}


swithListFun = function(paraList, ord) {
  res = list()
  for (i in 1:length(ord)) {
    res[[i]] = paraList[[ord[i]]]
  }
  return(res)

}

switchTaoFun = function(taoVec, ord) {
  res = c()
  for (i in 1:length(ord)) {
    res[i] = taoVec[ord[i]]
  }
  return(res)
}


switchZFun = function(label, ord) {
  res = c()
  for (i in 1:length(ord)) {
    res[label == i] = ord[i]
  }
  return(res)
}
