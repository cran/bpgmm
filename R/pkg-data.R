# set.seed(111)
#
# source("./R/sourceR.R")
#
#
# n <- 500 # sample size
# p <- 10 # number of variables
# q <- 4 # number of factors
# K <- 10 # true number of clusters
# # sINV_diag =  # diagonal of inverse variance of errors
#
# syntheticDataset <- simData(
#   sameLambda = TRUE, sameSigma = TRUE, K.true = K, n = n, q = q, p = p,
#   sINV_values = 1 / ((1:p))
# )
#
# ##
# nsim <- 10
# burn <- 20
#
# X <- t(syntheticDataset$data)
#
# # delta = 2; ggamma = 2
# # d_vec = c(1,1,1)
# # s_vec = c(1,1,1)
#
# q_new <- 4
#
# m_step <- 1
# v_step <- 1
# constraint <- c(0, 0, 0)
#
# m_init <- 20
# m_range <- c(1, 20)
#
# res <- pgmm_rjmcmc(X, m_init, m_range, q_new,
#   niter = nsim, burn = burn, constraint = constraint,
#   m_step = m_step, v_step = v_step
# )
#
# summarize_pgmm_rjmcmc(res, true_cluster = syntheticDataset$class)
