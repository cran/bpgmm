#' CalculateProposalLambda
#'
#' @param hparam hparam
#' @param thetaYList thetaYList
#' @param CxyList CxyList
#' @param constraint constraint
#' @param m the number of clusters
#' @param qVec the vector of the number of factors in each clusters
#' @param p the number of features
#'
#' @return calculated lambda for proposal function
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
#' constraint = c(1,1,1)
#'
#' hparam = new("Hparam", alpha1=1, alpha2=2, bbeta=3, delta=4, ggamma=5)
#' muBar = apply(X, MARGIN = 1, FUN = mean)
#' ZOneDim = kmeans(x = t(X), centers = m)$cluster
#' thetaYList = generatePriorThetaY(m, n, p, muBar, hparam, qVec, ZOneDim, constraint)
#'
#' \donttest{
#' url <- paste0("https://github.com/lzyacht/bpgmm-examples/",
#' "blob/master/data/CalculateProposalLambda.RData?raw=true")
#' download.file(url, destfile= "CalculateProposalLambda.RData", mode = "wb")
#' load("CalculateProposalLambda.RData")
#' CalculateProposalLambda(hparam, thetaYList, CxyList, constraint, m, p, qVec)
#' }
#'
#' @export

CalculateProposalLambda = function(hparam, thetaYList, CxyList, constraint, m, p, qVec){

    alpha1 = hparam@alpha1
    alpha2 = hparam@alpha2

    M = thetaYList@M
    psy = thetaYList@psy
    ##
    Cxxk = CxyList$Cxxk; Cxyk = CxyList$Cxyk; Cyyk = CxyList$Cyyk;
    Cytytk = CxyList$Cytytk; Cxtytk = CxyList$Cxtytk; CxL1k = CxyList$CxL1k;
    Cxmyk = CxyList$Cxmyk; sumCxmyk = CxyList$sumCxmyk; sumCyyk = CxyList$sumCyyk
    A = CxyList$A; nVec = CxyList$nVec

    ##
    lambda = list()
    if(constraint[1] == T & constraint[2] == T & constraint[3] == T){
        ##model 1

        sumCxmyk = 0
        sumCyyk = 0
        for(k in 1:m){
            sumCxmyk = sumCxmyk + Cxmyk[[k]]
            sumCyyk = sumCyyk + Cyyk[[k]] + alpha2 / m *diag(qVec[k])
        }

        for(k in 1:m){
            if(k == 1){
                lambda[[k]] = rmvnorm(1, mean = c(sumCxmyk %*% solve(sumCyyk)),
                sigma = kronecker(solve(sumCyyk), psy[[k]])
                )
                lambda[[k]] = matrix(lambda[[k]], p, qVec[k])

            }else{
                lambda[[k]] = lambda[[1]]
            }
        }
        ##model 1 end
    }else if(constraint[1] == T & constraint[2] == T & constraint[3] == F){
        ##model 2

        sumCxmyk = 0
        sumCyyk = 0
        for(k in 1:m){
            sumCxmyk = sumCxmyk + Cxmyk[[k]]
            sumCyyk = sumCyyk + Cyyk[[k]] + alpha2 / m *diag(qVec[k])
        }
        for(k in 1:m){
            if(k == 1){
                lambda[[k]] = rmvnorm(1, mean = c(sumCxmyk %*% solve(sumCyyk)),
                sigma =   kronecker(solve(sumCyyk), psy[[k]])
                )
                lambda[[k]] = matrix(lambda[[k]], p, qVec[k])

            }else{
                lambda[[k]] = lambda[[1]]
            }
        }

        ##end model 2
    }else if(constraint[1] == T & constraint[2] == F & constraint[3] == T){
        ##model 3


        sumPhiCxy = 0
        sumPhiCyy = 0

        for(k in 1:m){

            sumPhiCxy = sumPhiCxy + 1/psy[[k]][1,1] * Cxmyk[[k]]
            sumPhiCyy = sumPhiCyy + 1/psy[[k]][1,1] * (Cyyk[[k]] + alpha2 / m * diag(qVec[k]))
        }

        for(k in 1:m){
            if(k == 1){
                lambda[[k]] =   rmvnorm(1, mean = c(sumPhiCxy %*% solve(sumPhiCyy)),
                sigma = kronecker(solve(sumPhiCyy), diag(p)))

                lambda[[k]] = matrix(lambda[[k]], p, qVec[k])
            }else{
                lambda[[k]] = lambda[[1]]
            }
        }

        ## end model 3
    }else if(constraint[1] == T & constraint[2] == F & constraint[3] == F){
        ##model 4
        sumVar = 0
        B = 0
        for(k in 1:m){
            sumVar = sumVar + kronecker(Cyyk[[k]] + alpha2 /m * diag(qVec[k])
            , solve(psy[[k]]))
            B = B +  solve(psy[[k]]) %*% Cxmyk[[k]]

        }
        lambdaVar = solve(sumVar)
        lambdaMean = t(c(B)) %*% lambdaVar
        for(k in 1:m){
            if(k == 1){
                lambda[[k]] =   rmvnorm(1, mean = lambdaMean,
                sigma = lambdaVar)

                lambda[[k]] = matrix(lambda[[k]], p, qVec[k])
            }else{
                lambda[[k]] = lambda[[1]]
            }
        }
        ##end model 4
    }else if(constraint[1] == F & constraint[2] == T & constraint[3] == T){
        ##model 5
        for(k in 1:m){

            lambda[[k]] = rmvnorm(1, mean = c(Cxmyk[[k]] %*% solve(Cyyk[[k]]+  alpha2 * diag(qVec[k]))),
            sigma = kronecker(solve(Cyyk[[k]]+  alpha2 * diag(qVec[k])), psy[[k]])
            )
            lambda[[k]] = matrix(lambda[[k]], p, qVec[k])
        }

        ##end model 5
    }else if(constraint[1] == F & constraint[2] == T & constraint[3] == F){
        ##model 6
        for(k in 1:m){

            lambda[[k]] = rmvnorm(1, mean = c(Cxmyk[[k]] %*% solve(Cyyk[[k]]+  alpha2 * diag(qVec[k]))),
            sigma =  kronecker(solve(sumCyyk), psy[[k]])
            )
            lambda[[k]] = matrix(lambda[[k]], p, qVec[k])
        }

        ## end model 6
    }else if(constraint[1] == F & constraint[2] == F & constraint[3] == T){
        ##model 7
        for(k in 1:m){

            lambda[[k]] = rmvnorm(1, mean = c(Cxmyk[[k]] %*% solve(Cyyk[[k]]+  alpha2 * diag(qVec[k]))),
            sigma = kronecker(solve(sumCyyk), psy[[k]])
            )
            lambda[[k]] = matrix(lambda[[k]], p, qVec[k])
        }
        ## end model 7
    }else if(constraint[1] == F & constraint[2] == F & constraint[3] == F){
        ##model 8
        for(k in 1:m){

            lambda[[k]] = rmvnorm(1, mean = c(Cxmyk[[k]] %*% solve(Cyyk[[k]]+  alpha2 * diag(qVec[k]))),
            sigma = kronecker(solve(sumCyyk), psy[[k]])
            )
            lambda[[k]] = matrix(lambda[[k]], p, qVec[k])
        }
        ##end model 8
    }
    return(lambda)
}

#' EvaluateProposalLambda
#'
#' @param hparam hparam
#' @param thetaYList thetaYList
#' @param CxyList CxyList
#' @param constraint constraint
#' @param newlambda newlambda
#' @param m the number of clusters
#' @param qVec the vector of the number of factors in each clusters
#' @param p the number of features
#'
#' @return density value from evaluating lambda
#'
#' @export
#'


EvaluateProposalLambda = function(hparam, thetaYList, CxyList, constraint, newlambda, m, qVec, p){

    alpha1 = hparam@alpha1
    alpha2 = hparam@alpha2

    M = thetaYList@M
    psy = thetaYList@psy
    lambda = newlambda
    ##
    Cxxk = CxyList$Cxxk; Cxyk = CxyList$Cxyk; Cyyk = CxyList$Cyyk;
    Cytytk = CxyList$Cytytk; Cxtytk = CxyList$Cxtytk; CxL1k = CxyList$CxL1k;
    Cxmyk = CxyList$Cxmyk; sumCxmyk = CxyList$sumCxmyk; sumCyyk = CxyList$sumCyyk
    A = CxyList$A; nVec = CxyList$nVec

    ##
    lambdaEval = c()
    if(constraint[1] == T & constraint[2] == T & constraint[3] == T){
        ##model 1

        sumCxmyk = 0
        sumCyyk = 0
        for(k in 1:m){
            sumCxmyk = sumCxmyk + Cxmyk[[k]]
            sumCyyk = sumCyyk + Cyyk[[k]] + alpha2 / m *diag(qVec[k])
        }

        for(k in 1:m){
            if(k == 1){
                lambdaEval[k] = dmvnorm(x = c(lambda[[k]]),mean = c(sumCxmyk %*% solve(sumCyyk)),
                sigma = kronecker(solve(sumCyyk), psy[[k]]), log = T)
            }else{
                lambdaEval[k] =  0
            }
        }
        ##model 1 end
    }else if(constraint[1] == T & constraint[2] == T & constraint[3] == F){
        ##model 2

        sumCxmyk = 0
        sumCyyk = 0
        for(k in 1:m){
            sumCxmyk = sumCxmyk + Cxmyk[[k]]
            sumCyyk = sumCyyk + Cyyk[[k]] + alpha2 / m *diag(qVec[k])
        }

        for(k in 1:m){
            if(k == 1){
                lambdaEval[k] = dmvnorm(x = c(lambda[[k]]),mean = c(sumCxmyk %*% solve(sumCyyk)),
                sigma = kronecker(solve(sumCyyk), psy[[k]]),
                log = T)
            }else{
                lambdaEval[k] = 0
            }
        }

        ##end model 2
    }else if(constraint[1] == T & constraint[2] == F & constraint[3] == T){
        ##model 3

        sumPhiCxy = 0
        sumPhiCyy = 0

        for(k in 1:m){

            sumPhiCxy = sumPhiCxy + 1/psy[[k]][1,1] * Cxmyk[[k]]
            sumPhiCyy = sumPhiCyy + 1/psy[[k]][1,1] * (Cyyk[[k]] + alpha2 / m * diag(qVec[k]))
        }

        for(k in 1:m){
            if(k == 1){
                lambdaEval[k] =  dmvnorm(x = c(lambda[[k]]), mean = c(sumPhiCxy %*% solve(sumPhiCyy)),
                sigma = kronecker(solve(sumPhiCyy), diag(p)),log = T)


            }else{
                lambdaEval[k] = 0
            }
        }
        ## end model 3
    }else if(constraint[1] == T & constraint[2] == F & constraint[3] == F){
        ##model 4
        sumVar = 0
        B = 0
        for(k in 1:m){
            sumVar = sumVar + kronecker(Cyyk[[k]] + alpha2 /m * diag(qVec[k])
            , solve(psy[[k]]))
            B = B +  solve(psy[[k]]) %*% Cxmyk[[k]]

        }
        lambdaVar = solve(sumVar)
        lambdaMean = t(c(B)) %*% lambdaVar
        for(k in 1:m){
            if(k == 1){
                lambdaEval[k] = dmvnorm(x = c(lambda[[k]]), mean = lambdaMean,
                sigma = lambdaVar, log = T)
            }else{
                lambdaEval[k] = 0
            }
        }
        ##end model 4
    }else if(constraint[1] == F & constraint[2] == T & constraint[3] == T){
        ##model 5
        for(k in 1:m){
            lambdaEval[k] = dmvnorm(x = c(lambda[[k]]), mean = c(Cxmyk[[k]] %*% solve(Cyyk[[k]]+  alpha2 * diag(qVec[k]))),
            sigma = kronecker(solve(sumCyyk), psy[[k]]),log = T)
        }

        ##end model 5
    }else if(constraint[1] == F & constraint[2] == T & constraint[3] == F){
        ##model 6
        for(k in 1:m){
            lambdaEval[k] = dmvnorm(x = c(lambda[[k]]), mean = c(Cxmyk[[k]] %*% solve(Cyyk[[k]]+  alpha2 * diag(qVec[k]))),
            sigma = kronecker(solve(sumCyyk), psy[[k]]), log = T)
        }

        ## end model 6
    }else if(constraint[1] == F & constraint[2] == F & constraint[3] == T){
        ##model 7
        for(k in 1:m){

            lambdaEval[k] = dmvnorm(x = c(lambda[[k]]), mean = c(Cxmyk[[k]] %*% solve(Cyyk[[k]]+  alpha2 * diag(qVec[k]))),
            sigma = kronecker(solve(sumCyyk), psy[[k]]), log = T)
        }
        ## end model 7
    }else if(constraint[1] == F & constraint[2] == F & constraint[3] == F){
        ##model 8
        for(k in 1:m){

            lambdaEval[k] = dmvnorm(c(lambda[[k]]), mean = c(Cxmyk[[k]] %*% solve(Cyyk[[k]]+ alpha2 * diag(qVec[k]))),
            sigma = kronecker(solve(sumCyyk), psy[[k]]), log = T)
        }
        ##end model 8
    }
    return(sum(lambdaEval))
}


