#' Summarize RJMCMC Samples from a Bayesian PGMM Fit
#'
#' Summarizes posterior samples from [pgmm_rjmcmc()] into the modal allocation,
#' posterior counts for the number of clusters, posterior counts for the eight
#' PGMM covariance-constraint models, and optionally the adjusted Rand index
#' against a known reference partition.
#'
#' @param fit Result list from [pgmm_rjmcmc()].
#' @param true_cluster Optional true or reference cluster allocation.
#'
#' @return A list with `allocation`, `n_clusters`, `n_constraints`, and optionally
#'   `ari`.
#' @importFrom mclust adjustedRandIndex
#' @name summarize_pgmm_rjmcmc
#' @export
summarize_pgmm_rjmcmc <- function(fit, true_cluster = NULL) {
  if (!is.list(fit)) {
    stop("fit must be a result list from pgmm_rjmcmc()", call. = FALSE)
  }
  if (!is.list(fit$allocation_samples) || length(fit$allocation_samples) == 0L) {
    stop("fit$allocation_samples must contain at least one allocation sample", call. = FALSE)
  }
  if (!is.list(fit$constraint_samples) || length(fit$constraint_samples) != length(fit$allocation_samples)) {
    stop("fit$constraint_samples must match fit$allocation_samples", call. = FALSE)
  }

  allocation <- summarize_allocations(fit$allocation_samples)

  n_clusters <- table(sapply(fit$allocation_samples, function(x) {
    length(unique(x))
  }))

  n_constraints <- fit$constraint_samples
  n_constraints <- constraint_list_to_models(n_constraints)
  n_constraints <- table(n_constraints, dnn = "")



  summary <- list(allocation = allocation, n_clusters = n_clusters, n_constraints = n_constraints)

  if (!is.null(true_cluster)) {
    ari <- adjustedRandIndex(true_cluster, allocation)
    summary$ari <- ari
  }

  summary
}
