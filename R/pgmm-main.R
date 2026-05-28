#' Bayesian Model-Based Clustering with Parsimonious Gaussian Mixture Models
#'
#' Carries out model-based clustering using parsimonious Gaussian mixture
#' models. MCMC is used for parameter estimation and RJMCMC is used for model
#' selection.
#'
#' The `constraint` argument follows the three-letter PGMM model notation used
#' in Lu, Li, and Love (2021). The first entry indicates whether loading
#' matrices are shared across clusters, the second whether noise covariance
#' matrices are shared across clusters, and the third whether the noise
#' covariance is isotropic within each cluster. Use [model_to_constraint()] to
#' convert model names such as `CCC`, `CCU`, `CUC`, `CUU`, `UCC`, `UCU`, `UUC`,
#' and `UUU` into the numeric vector used internally.
#'
#' @import stats MASS mcmcse pgmm label.switching fabMix
#' @param X the observation matrix with variables in rows and observations in
#'   columns.
#' @param m_init the number of initial clusters.
#' @param m_range the allowed range for the number of clusters.
#' @param q_new the number of latent factors for a new cluster.
#' @param delta scalar hyperparameter for the noise covariance prior
#' @param ggamma scalar hyperparameter used in covariance-structure proposals
#' @param burn the number of burn-in iterations
#' @param niter the number of posterior sampling iterations
#' @param constraint initial PGMM covariance constraint. Use a three-letter
#'   model label such as `"CCC"` or `"UUU"`, or a numeric vector of length
#'   three with binary entries. For example, `c(1, 1, 1)` is `CCC`, the fully
#'   constrained model, and `c(0, 0, 0)` is `UUU`, the fully unconstrained
#'   model.
#' @param d_vec a vector of hyperparameters with length three, shape parameters
#'   for alpha1, alpha2 and bbeta respectively
#' @param s_vec a vector of hyperparameters with length three, rate parameters
#'   for alpha1, alpha2 and bbeta respectively
#' @param m_step indicator for RJMCMC model selection on the number of clusters.
#' @param v_step indicator for RJMCMC model selection on covariance structures.
#' @param split_combine indicator for using split/combine moves in the
#'   cluster-number RJMCMC step.
#' @param verbose logical; if `TRUE`, print iteration progress.
#'
#' @return A list of posterior samples with snake_case fields:
#'   `tau_samples`, `psi_samples`, `mean_samples`, `lambda_samples`,
#'   `factor_score_samples`, `allocation_samples`, `constraint_samples`,
#'   `alpha1_samples`, `alpha2_samples`, `beta_samples`, and
#'   `active_cluster_samples`.
#' @name pgmm_rjmcmc
#' @export

pgmm_rjmcmc <- function(X,
                       m_init,
                       m_range,
                       q_new,
                       delta = 2,
                       ggamma = 2,
                       burn = 20,
                       niter = 1000,
                       constraint = c(0, 0, 0),
                       d_vec = c(1, 1, 1),
                       s_vec = c(1, 1, 1),
                       m_step = 0,
                       v_step = 0,
                       split_combine = 0,
                       verbose = TRUE) {
  if (!is.matrix(X) || !is.numeric(X)) {
    stop("X must be a numeric matrix with variables in rows and observations in columns", call. = FALSE)
  }
  if (any(!is.finite(X))) {
    stop("X must contain only finite numeric values", call. = FALSE)
  }
  if (!isTRUE(length(m_range) == 2L) || any(!is.finite(m_range)) || any(m_range < 1) || m_range[1] > m_range[2]) {
    stop("m_range must be a length-two increasing positive numeric vector", call. = FALSE)
  }
  if (!isTRUE(length(m_init) == 1L) || !is.finite(m_init) || m_init < 1) {
    stop("m_init must be a positive scalar", call. = FALSE)
  }
  if (m_init < m_range[1] || m_init > m_range[2]) {
    stop("m_init must lie within m_range", call. = FALSE)
  }
  if (m_init > ncol(X)) {
    stop("m_init cannot exceed the number of observations in X", call. = FALSE)
  }
  if (!isTRUE(length(q_new) == 1L) || !is.finite(q_new) || q_new < 1) {
    stop("q_new must be a positive scalar", call. = FALSE)
  }
  if (!isTRUE(length(burn) == 1L) || !is.finite(burn) || burn < 0) {
    stop("burn must be a non-negative scalar", call. = FALSE)
  }
  if (!isTRUE(length(niter) == 1L) || !is.finite(niter) || niter < 0) {
    stop("niter must be a non-negative scalar", call. = FALSE)
  }
  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
    stop("verbose must be TRUE or FALSE", call. = FALSE)
  }
  if (!isTRUE(length(m_step) == 1L) || is.na(m_step) || !(m_step %in% c(0, 1))) {
    stop("m_step must be 0 or 1", call. = FALSE)
  }
  if (!isTRUE(length(v_step) == 1L) || is.na(v_step) || !(v_step %in% c(0, 1))) {
    stop("v_step must be 0 or 1", call. = FALSE)
  }
  if (!isTRUE(length(split_combine) == 1L) || is.na(split_combine) || !(split_combine %in% c(0, 1))) {
    stop("split_combine must be 0 or 1", call. = FALSE)
  }

  m_init <- as.integer(m_init)
  m_range <- as.integer(m_range)
  q_new <- as.integer(q_new)
  burn <- as.integer(burn)
  niter <- as.integer(niter)
  m_step <- as.integer(m_step)
  v_step <- as.integer(v_step)
  split_combine <- as.integer(split_combine)
  if (is.character(constraint)) {
    constraint <- model_to_constraint(constraint)
  } else {
    constraint <- model_to_constraint(constraint_to_model(constraint))
  }

  n <- ncol(X)
  p <- nrow(X)

  alpha1 <- rgamma(1, d_vec[1], s_vec[1])
  alpha2 <- rgamma(1, d_vec[2], s_vec[2])
  bbeta <- rgamma(1, d_vec[3], s_vec[3])

  hparam <- new("Hparam", alpha1 = alpha1, alpha2 = alpha2, bbeta = bbeta, delta = delta, ggamma = ggamma)

  hparamInit <- hparam

  muBar <- X[, sample.int(n, 1)]

  ## cluster indicator
  clusInd <- rep(0, m_range[2])
  clusInd[1:m_init] <- 1

  ## qinit
  qVec <- rep(0, m_range[2])
  qVec[1:m_init] <- q_new


  ## priors
  ZOneDim <- kmeans(x = t(X), centers = m_init)$cluster
  thetaYList <- generate_prior_theta_y(m_init, n, p, muBar, hparam, qVec, ZOneDim, constraint)

  ## burn in
  for (i in seq_len(burn)) {
    MCMCobj <- stay_mcmc_update(X, thetaYList, ZOneDim, hparam, qVec, q_new, d_vec, s_vec, constraint, clusInd)
    ZOneDim <- MCMCobj$ZOneDim
    thetaYList <- MCMCobj$thetaYList
    hparam <- MCMCobj$hparam
    hparam@alpha2 <- max(0.01, hparam@alpha2)
  }

  thetaYList <- clear_inactive_components(thetaYList, clusInd, m_range[2])
  ##
  alpha1_samples <- c()
  alpha2_samples <- c()
  beta_samples <- c()
  tau_samples <- list()
  psi_samples <- list()
  mean_samples <- list()
  lambda_samples <- list()
  factor_score_samples <- list()
  allocation_samples <- list()
  constraint_samples <- list()
  active_cluster_samples <- list()
  ##

  for (h in seq_len(niter)) {
    if (verbose) {
      cat("iter = ", h, "======>\n")
    }

    ## choose m or choose v
    if (m_step == 1) {
      MCMCobj <- mstep_rjmcmc_update(X, muBar, p, thetaYList, ZOneDim, hparam, hparamInit, qVec, q_new, d_vec, s_vec, constraint, clusInd, m_range, "BD")
      ZOneDim <- MCMCobj$ZOneDim
      thetaYList <- MCMCobj$thetaYList
      hparam <- MCMCobj$hparam
      qVec <- MCMCobj$qVec
      clusInd <- MCMCobj$clusInd
      ##
      if (split_combine == 1) {
        MCMCobj <- mstep_rjmcmc_update(X, muBar, p, thetaYList, ZOneDim, hparam, hparamInit, qVec, q_new, d_vec, s_vec, constraint, clusInd, m_range, "SC")
        ZOneDim <- MCMCobj$ZOneDim
        thetaYList <- MCMCobj$thetaYList
        hparam <- MCMCobj$hparam
        qVec <- MCMCobj$qVec
        clusInd <- MCMCobj$clusInd
      }
    }

    if (v_step == 1) {
      MCMCobj <- vstep_rjmcmc_update(X, muBar, p, thetaYList, ZOneDim, hparam, hparamInit, qVec, q_new, ggamma, d_vec, s_vec, constraint, clusInd)
      ZOneDim <- MCMCobj$ZOneDim
      thetaYList <- MCMCobj$thetaYList
      hparam <- MCMCobj$hparam
      constraint <- MCMCobj$constraint
    }

    # stay step
    MCMCobj <- stay_mcmc_update(X, thetaYList, ZOneDim, hparam, qVec, q_new, d_vec, s_vec, constraint, clusInd)
    ZOneDim <- MCMCobj$ZOneDim
    thetaYList <- MCMCobj$thetaYList
    hparam <- MCMCobj$hparam
    hparam@alpha2 <- max(0.01, hparam@alpha2)

    ## save
    active_cluster_samples[[h]] <- clusInd
    alpha1_samples[h] <- hparam@alpha1
    alpha2_samples[h] <- hparam@alpha2
    beta_samples[h] <- hparam@bbeta
    tau_samples[[h]] <- thetaYList@tao
    psi_samples[[h]] <- thetaYList@psy
    mean_samples[[h]] <- thetaYList@M
    lambda_samples[[h]] <- thetaYList@lambda
    factor_score_samples[[h]] <- thetaYList@Y
    allocation_samples[[h]] <- ZOneDim
    constraint_samples[[h]] <- constraint
  }

  list(
    tau_samples = tau_samples,
    psi_samples = psi_samples,
    mean_samples = mean_samples,
    lambda_samples = lambda_samples,
    factor_score_samples = factor_score_samples,
    allocation_samples = allocation_samples,
    constraint_samples = constraint_samples,
    alpha1_samples = alpha1_samples,
    alpha2_samples = alpha2_samples,
    beta_samples = beta_samples,
    active_cluster_samples = active_cluster_samples
  )
}
