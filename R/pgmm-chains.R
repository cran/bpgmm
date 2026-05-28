#' Run Multiple Independent Bayesian PGMM Chains
#'
#' Runs independent `pgmm_rjmcmc()` chains, optionally in parallel. This is the
#' safest way to use multiple CPU cores because each MCMC iteration depends on
#' the previous state, while independent chains can be evaluated concurrently.
#'
#' @inheritParams pgmm_rjmcmc
#' @param chains positive integer giving the number of independent chains.
#' @param cores positive integer giving the number of worker processes to use.
#'   Values greater than `chains` are reduced to `chains`.
#' @param seed optional integer seed used to generate deterministic per-chain
#'   seeds.
#'
#' @return A list with one fitted `pgmm_rjmcmc()` result per chain. The result
#'   has class `bpgmm_rjmcmc_chains` and stores the per-chain seeds in the
#'   `chain_seeds` attribute.
#' @export
pgmm_rjmcmc_chains <- function(X,
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
                               verbose = FALSE,
                               chains = 2,
                               cores = min(chains, available_cores()),
                               seed = NULL) {
  chains <- validate_positive_integer(chains, "chains")
  cores <- validate_positive_integer(cores, "cores")
  cores <- min(cores, chains)

  if (!is.null(seed)) {
    if (!isTRUE(length(seed) == 1L) || !is.finite(seed)) {
      stop("seed must be NULL or a finite scalar", call. = FALSE)
    }
    old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    } else {
      NULL
    }
    on.exit(restore_random_seed(old_seed), add = TRUE)
    set.seed(as.integer(seed))
  }

  chain_seeds <- sample.int(.Machine$integer.max, chains)
  run_chain <- function(chain_id) {
    set.seed(chain_seeds[[chain_id]])
    pgmm_rjmcmc(
      X = X,
      m_init = m_init,
      m_range = m_range,
      q_new = q_new,
      delta = delta,
      ggamma = ggamma,
      burn = burn,
      niter = niter,
      constraint = constraint,
      d_vec = d_vec,
      s_vec = s_vec,
      m_step = m_step,
      v_step = v_step,
      split_combine = split_combine,
      verbose = verbose
    )
  }

  fits <- if (cores == 1L) {
    lapply(seq_len(chains), run_chain)
  } else if (.Platform$OS.type == "windows") {
    cluster <- parallel::makeCluster(cores)
    on.exit(parallel::stopCluster(cluster), add = TRUE)
    parallel::clusterEvalQ(cluster, library(bpgmm))
    parallel::clusterExport(
      cluster,
      varlist = c(
        "X", "m_init", "m_range", "q_new", "delta", "ggamma", "burn",
        "niter", "constraint", "d_vec", "s_vec", "m_step", "v_step",
        "split_combine", "verbose", "chain_seeds"
      ),
      envir = environment()
    )
    parallel::parLapply(cluster, seq_len(chains), function(chain_id) {
      set.seed(chain_seeds[[chain_id]])
      bpgmm::pgmm_rjmcmc(
        X = X,
        m_init = m_init,
        m_range = m_range,
        q_new = q_new,
        delta = delta,
        ggamma = ggamma,
        burn = burn,
        niter = niter,
        constraint = constraint,
        d_vec = d_vec,
        s_vec = s_vec,
        m_step = m_step,
        v_step = v_step,
        split_combine = split_combine,
        verbose = verbose
      )
    })
  } else {
    parallel::mclapply(seq_len(chains), run_chain, mc.cores = cores)
  }

  names(fits) <- paste0("chain_", seq_len(chains))
  attr(fits, "chain_seeds") <- chain_seeds
  class(fits) <- c("bpgmm_rjmcmc_chains", class(fits))
  fits
}
