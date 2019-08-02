#' CalculateProposalPsy
#'
#' @param hparam hparam
#' @param thetaYList thetaYList
#' @param CxyList CxyList
#' @param constraint constraint
#' @param m the number of clusters
#' @param qVec the vector of the number of factors in each clusters
#' @param p the number of features
#'
#' @return calculated psy for proposal function
#'
#' @export

CalculateProposalPsy = function(hparam, thetaYList, CxyList, constraint, m, p, qVec){

    alpha1 = hparam@alpha1
    alpha2 = hparam@alpha2
    bbeta  = hparam@bbeta
    delta  = hparam@delta
    M      = thetaYList@M
    psy    = thetaYList@psy
    lambda = thetaYList@lambda


    Cxxk = CxyList$Cxxk; Cxyk = CxyList$Cxyk; Cyyk = CxyList$Cyyk;
    Cytytk = CxyList$Cytytk; Cxtytk = CxyList$Cxtytk; CxL1k = CxyList$CxL1k;
    Cxmyk = CxyList$Cxmyk; sumCxmyk = CxyList$sumCxmyk; sumCyyk = CxyList$sumCyyk
    A = CxyList$A; nVec = CxyList$nVec

    # post tilda lambda_k = {mu_k, lambda_k}, first column is mu_k
    tildaLambda = list()
    for(k in 1:m){
        tildaLambda[[k]] = cbind(t(M[[k]]), lambda[[k]])
    }

    if(constraint[1] == T & constraint[2] == T & constraint[3] == T){
        ##model 1

        ## post psy
        psy = list()

        shapePara = 0
        ratePara = 0
        for(k in 1:m){
            shapePara = shapePara + p/2 * (nVec[k] + qVec[k]/m + (2 * delta - 2)/(m*p)  + 1)
            ratePara = ratePara + 1/2 * diag(Cxxk[[k]] - 2 * Cxtytk[[k]] %*% t(tildaLambda[[k]])
            + tildaLambda[[k]]%*%(Cytytk[[k]] + A[[k]]/m)%*%t(tildaLambda[[k]] )
            + (2 * bbeta)/(m*p) * diag(rep(1,p)))

        }
        shapePara = shapePara + 1
        ratePara = sum(ratePara)

        invpsy = rgamma(1, shape = shapePara, rate = ratePara)
        for(k in 1:m){
            psy[[k]] = diag(rep(1/invpsy, p))
        }
        ##model 1 end
    }else if(constraint[1] == T & constraint[2] == T & constraint[3] == F){
        ##model 2

        ## post psy
        psy = list()

        shapePara = 0
        ratePara = 0
        for(k in 1:m){
            shapePara = shapePara + 1/2 * (nVec[k] + qVec[k]/m + (2 * delta - 2)/m + 1)
            ratePara = ratePara + 1/2 * diag(Cxxk[[k]] - 2 * Cxtytk[[k]] %*% t(tildaLambda[[k]])
            + tildaLambda[[k]]%*%(Cytytk[[k]] + A[[k]]/m)%*%t(tildaLambda[[k]] )
            + 2 * bbeta/m * diag(rep(1,p)))

        }
        shapePara = shapePara + 1

        invpsy = c()
        for(j in 1:p){
            invpsy[j] = rgamma(1, shape = shapePara, rate = ratePara[j])
        }
        for(k in 1:m){
            psy[[k]] = diag(1/invpsy)
        }
        ##end model 2
    }else if(constraint[1] == T & constraint[2] == F & constraint[3] == T){
        ##model 3

        ## post psy
        psy = list()

        shapePara = 0
        ratePara = 0
        for(k in 1:m){
            shapePara = p/2 * (nVec[k] + qVec[k]/m + (2 * delta - 2)/p + 1) + 1
            ratePara = 1/2 * sum(diag(Cxxk[[k]] - 2 * Cxtytk[[k]] %*% t(tildaLambda[[k]])
            + tildaLambda[[k]]%*%(Cytytk[[k]] + A[[k]])%*%t(tildaLambda[[k]] )
            + 2 * bbeta /p  * diag(rep(1,p))))

            invpsy = rgamma(1, shape = shapePara, rate = ratePara)
            psy[[k]] = diag(rep(1/invpsy,p))
        }


        ## end model 3
    }else if(constraint[1] == T & constraint[2] == F & constraint[3] == F){
        ##model 4
        ## post psy
        psy = list()

        shapePara = 0
        ratePara = 0
        for(k in 1:m){
            shapePara = 1/2 * (nVec[k] + qVec[k]/m + 2 * delta - 1) + 1
            ratePara = 1/2 * diag(Cxxk[[k]] - 2 * Cxtytk[[k]] %*% t(tildaLambda[[k]])
            + tildaLambda[[k]]%*%(Cytytk[[k]] + A[[k]]/m)%*%t(tildaLambda[[k]] )
            + 2 * bbeta * diag(rep(1,p)))


            invpsy = c()
            for(j in 1:p){
                invpsy[j] = rgamma(1, shape = shapePara, rate = ratePara[j])
            }
            # invpsy = rgamma(p, shape = shapePara, rate = ratePara)
            psy[[k]] = diag(1/invpsy)
        }

        ##end model 4
    }else if(constraint[1] == F & constraint[2] == T & constraint[3] == T){
        ##model 5
        ## post psy
        psy = list()

        shapePara = 0
        ratePara = 0
        for(k in 1:m){
            shapePara = shapePara + p/2 * (nVec[k] + qVec[k] + (2 * delta - 2)/(m*p) + 1)
            ratePara = ratePara + 1/2 * diag(Cxxk[[k]] - 2 * Cxtytk[[k]] %*% t(tildaLambda[[k]])
            + tildaLambda[[k]]%*%(Cytytk[[k]] + A[[k]])%*%t(tildaLambda[[k]] )
            + 2 * bbeta/(m*p) * diag(rep(1,p)))

        }
        shapePara = shapePara + 1
        ratePara = sum(ratePara)

        invpsy = rgamma(1, shape = shapePara, rate = ratePara)

        for(k in 1:m){
            psy[[k]] = diag(rep(1/invpsy),p)
        }
        ##end model 5
    }else if(constraint[1] == F & constraint[2] == T & constraint[3] == F){
        ##model 6
        ## post psy
        psy = list()

        shapePara = 0
        ratePara = 0
        for(k in 1:m){
            shapePara = shapePara + 1/2 * (nVec[k] + qVec[k] + (2 * delta - 2)/m + 1)
            ratePara = ratePara + 1/2 * diag(Cxxk[[k]] - 2 * Cxtytk[[k]] %*% t(tildaLambda[[k]])
            + tildaLambda[[k]]%*%(Cytytk[[k]] + A[[k]])%*%t(tildaLambda[[k]] )
            + 2 * bbeta /m  * diag(rep(1,p)))

        }
        shapePara = shapePara + 1

        invpsy = c()
        for(j in 1:p){
            invpsy[j] = rgamma(1, shape = shapePara, rate = ratePara[j])
        }
        for(k in 1:m){
            psy[[k]] = diag(1/invpsy)
        }
        ## end model 6
    }else if(constraint[1] == F & constraint[2] == F & constraint[3] == T){
        ##model 7

        ## post psy
        psy = list()

        shapePara = 0
        ratePara = 0
        for(k in 1:m){
            shapePara = p/2 * (nVec[k] + qVec[k] + (2 * delta - 2)/p + 1) + 1
            ratePara = 1/2 * sum(diag(Cxxk[[k]] - 2 * Cxtytk[[k]] %*% t(tildaLambda[[k]])
            + tildaLambda[[k]]%*%(Cytytk[[k]] + A[[k]])%*%t(tildaLambda[[k]] )
            + 2 * bbeta /p * diag(rep(1,p))))

            invpsy = rgamma(1, shape = shapePara, rate = ratePara)
            psy[[k]] = diag(rep(1/invpsy,p))
        }
        ## end model 7
    }else if(constraint[1] == F & constraint[2] == F & constraint[3] == F){
        ##model 8
        ## post psy
        psy = list()

        shapePara = 0
        ratePara = 0
        for(k in 1:m){
            shapePara = 1/2 * (nVec[k] + qVec[k] + 2 * delta - 1) + 1
            ratePara = 1/2 * diag(Cxxk[[k]] - 2 * Cxtytk[[k]] %*% t(tildaLambda[[k]])
            + tildaLambda[[k]]%*%(Cytytk[[k]] + A[[k]])%*%t(tildaLambda[[k]] )
            + 2 * bbeta * diag(rep(1,p)))

            invpsy = rgamma(p, shape = shapePara, rate = ratePara)
            psy[[k]] = diag(1/invpsy)
        }

        ##end model 8
    }
    return(psy)
}


#' EvaluateProposalPsy
#'
#' @param hparam hparam
#' @param thetaYList thetaYList
#' @param CxyList CxyList
#' @param constraint constraint
#' @param newpsy newpsy
#' @param m the number of clusters
#' @param qVec the vector of the number of factors in each clusters
#' @param p the number of features
#'
#' @return evaluation of proposal psy
#' @export
#'

EvaluateProposalPsy = function(hparam, thetaYList, CxyList, constraint, newpsy, m, p, qVec){

    alpha1 = hparam@alpha1
    alpha2 = hparam@alpha2
    bbeta = hparam@bbeta
    delta = hparam@delta
    M = thetaYList@M
    lambda = thetaYList@lambda

    ##
    Cxxk = CxyList$Cxxk; Cxyk = CxyList$Cxyk; Cyyk = CxyList$Cyyk;
    Cytytk = CxyList$Cytytk; Cxtytk = CxyList$Cxtytk; CxL1k = CxyList$CxL1k;
    Cxmyk = CxyList$Cxmyk; sumCxmyk = CxyList$sumCxmyk; sumCyyk = CxyList$sumCyyk
    A = CxyList$A; nVec = CxyList$nVec

    ##
    psyEval = c()
    psy = newpsy
    tildaLambda = list()
    for(k in 1:m){
        tildaLambda[[k]] = cbind(t(M[[k]]), lambda[[k]])
    }
    ##
    if(constraint[1] == T & constraint[2] == T & constraint[3] == T){
        ##model 1

        ## post psy
        shapePara = 0
        ratePara = 0
        for(k in 1:m){
            shapePara = shapePara + p/2 * (nVec[k] + qVec[k]/m + (2 * delta - 2)/(m*p)  + 1)
            ratePara = ratePara + 1/2 * diag(Cxxk[[k]] - 2 * Cxtytk[[k]] %*% t(tildaLambda[[k]])
            + tildaLambda[[k]]%*%(Cytytk[[k]] + A[[k]]/m)%*%t(tildaLambda[[k]] )
            + (2 * bbeta)/(m*p) * diag(rep(1,p)))

        }
        shapePara = shapePara + 1
        ratePara = sum(ratePara)

        psyEval = dgamma(1/psy[[1]][1,1], shape = shapePara, rate = ratePara, log = T)
        ##model 1 end
    }else if(constraint[1] == T & constraint[2] == T & constraint[3] == F){
        ##model 2

        shapePara = 0
        ratePara = 0
        for(k in 1:m){
            shapePara = shapePara + 1/2 * (nVec[k] + qVec[k]/m + (2 * delta - 2)/m + 1)
            ratePara = ratePara + 1/2 * diag(Cxxk[[k]] - 2 * Cxtytk[[k]] %*% t(tildaLambda[[k]])
            + tildaLambda[[k]]%*%(Cytytk[[k]] + A[[k]]/m)%*%t(tildaLambda[[k]] )
            + 2 * bbeta/m * diag(rep(1,p)))

        }
        shapePara = shapePara + 1

        invpsy = 1/diag(psy[[1]])
        for(j in 1:p){
            psyEval[j] = dgamma(invpsy[j], shape = shapePara, rate = ratePara[j], log = T)
        }
        ##end model 2
    }else if(constraint[1] == T & constraint[2] == F & constraint[3] == T){
        ## post psy

        shapePara = 0
        ratePara = 0
        for(k in 1:m){
            shapePara = p/2 * (nVec[k] + qVec[k]/m + (2 * delta - 2)/p + 1) + 1
            ratePara = 1/2 * sum(diag(Cxxk[[k]] - 2 * Cxtytk[[k]] %*% t(tildaLambda[[k]])
            + tildaLambda[[k]]%*%(Cytytk[[k]] + A[[k]])%*%t(tildaLambda[[k]] )
            + 2 * bbeta /p  * diag(rep(1,p))))

            psyEval[k] = dgamma(1/psy[[k]][1,1], shape = shapePara, rate = ratePara, log = T)
        }


        ## end model 3
    }else if(constraint[1] == T & constraint[2] == F & constraint[3] == F){
        ##model 4
        ## post psy

        shapePara = 0
        ratePara = 0
        psyEval = matrix(NA, nrow = m, ncol = p)
        for(k in 1:m){
            shapePara = 1/2 * (nVec[k] + qVec[k]/m + 2 * delta - 1) + 1
            ratePara = 1/2 * diag(Cxxk[[k]] - 2 * Cxtytk[[k]] %*% t(tildaLambda[[k]])
            + tildaLambda[[k]]%*%(Cytytk[[k]] + A[[k]]/m)%*%t(tildaLambda[[k]] )
            + 2 * bbeta * diag(rep(1,p)))


            invpsy = 1/diag(psy[[k]])
            for(j in 1:p){
                psyEval[k,j] = dgamma(invpsy[j], shape = shapePara, rate = ratePara[j], log = T)
            }
        }

        ##end model 4
    }else if(constraint[1] == F & constraint[2] == T & constraint[3] == T){
        ##model 5

        ## post psy
        shapePara = 0
        ratePara = 0
        for(k in 1:m){
            shapePara = shapePara + p/2 * (nVec[k] + qVec[k] + (2 * delta - 2)/(m*p) + 1)
            ratePara = ratePara + 1/2 * diag(Cxxk[[k]] - 2 * Cxtytk[[k]] %*% t(tildaLambda[[k]])
            + tildaLambda[[k]]%*%(Cytytk[[k]] + A[[k]])%*%t(tildaLambda[[k]] )
            + 2 * bbeta/(m*p) * diag(rep(1,p)))

        }
        shapePara = shapePara + 1
        ratePara = sum(ratePara)

        psyEval = dgamma(1/psy[[1]][1,1], shape = shapePara, rate = ratePara, log = T)

        ##end model 5
    }else if(constraint[1] == F & constraint[2] == T & constraint[3] == F){
        ##model 6


        ## post psy

        shapePara = 0
        ratePara = 0
        for(k in 1:m){
            shapePara = shapePara + 1/2 * (nVec[k] + qVec[k] + (2 * delta - 2)/m + 1)
            ratePara = ratePara + 1/2 * diag(Cxxk[[k]] - 2 * Cxtytk[[k]] %*% t(tildaLambda[[k]])
            + tildaLambda[[k]]%*%(Cytytk[[k]] + A[[k]])%*%t(tildaLambda[[k]] )
            + 2 * bbeta /m  * diag(rep(1,p)))

        }
        shapePara = shapePara + 1

        invpsy = 1/diag(psy[[1]])
        for(j in 1:p){
            psyEval[j] = dgamma(invpsy[j], shape = shapePara, rate = ratePara[j], log = T)
        }

        ## end model 6
    }else if(constraint[1] == F & constraint[2] == F & constraint[3] == T){
        ##model 7


        ## post psy

        shapePara = 0
        ratePara = 0
        for(k in 1:m){
            shapePara = p/2 * (nVec[k] + qVec[k] + (2 * delta - 2)/p + 1) + 1
            ratePara = 1/2 * sum(diag(Cxxk[[k]] - 2 * Cxtytk[[k]] %*% t(tildaLambda[[k]])
            + tildaLambda[[k]]%*%(Cytytk[[k]] + A[[k]])%*%t(tildaLambda[[k]] )
            + 2 * bbeta /p * diag(rep(1,p))))

            psyEval[k] = dgamma(1/psy[[k]][1,1], shape = shapePara, rate = ratePara, log = T)
        }
        ## end model 7
    }else if(constraint[1] == F & constraint[2] == F & constraint[3] == F){
        ##model 8

        ## post psy
        psyEval = matrix(NA, nrow = m, ncol = p)
        shapePara = 0
        ratePara = 0
        for(k in 1:m){
            shapePara = 1/2 * (nVec[k] + qVec[k] + 2 * delta - 1) + 1
            ratePara = 1/2 * diag(Cxxk[[k]] - 2 * Cxtytk[[k]] %*% t(tildaLambda[[k]])
            + tildaLambda[[k]]%*%(Cytytk[[k]] + A[[k]])%*%t(tildaLambda[[k]] )
            + 2 * bbeta * diag(rep(1,p)))


            invpsy = 1/diag(psy[[k]])
            for(j in 1:p){
                psyEval[k,j] = dgamma(invpsy[j], shape = shapePara, rate = ratePara[j], log = T)
            }
        }
        ##end model 8
    }
    return(sum(psyEval))
}
