## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library(bpgmm)

set.seed(2029)
X <- cbind(
  matrix(rnorm(10, mean = -2.0, sd = 0.25), nrow = 2),
  matrix(rnorm(10, mean =  0.0, sd = 0.25), nrow = 2),
  matrix(rnorm(10, mean =  2.0, sd = 0.25), nrow = 2)
)
known_labels <- rep(1:3, each = 5)

## ----fig.width = 5.5, fig.height = 4, fig.alt = "Scatter plot of a compact three-cluster diagnostic data set."----
cluster_cols <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7")
plot(
  X[1, ], X[2, ],
  col = cluster_cols[known_labels],
  pch = 19,
  xlab = "Variable 1",
  ylab = "Variable 2",
  main = "Diagnostic example",
  asp = 1
)

## -----------------------------------------------------------------------------
fits <- pgmm_rjmcmc_chains(
  X = X,
  m_init = 3,
  m_range = c(1, 4),
  q_new = 1,
  burn = 1,
  niter = 4,
  constraint = "UUU",
  m_step = 1,
  v_step = 1,
  chains = 2,
  cores = 1,
  seed = 2029,
  verbose = FALSE
)

length(fits)
attr(fits, "chain_seeds")

## -----------------------------------------------------------------------------
chain_summaries <- lapply(fits, summarize_pgmm_rjmcmc, true_cluster = known_labels)

data.frame(
  chain = names(chain_summaries),
  ari = vapply(chain_summaries, function(x) x$ari, numeric(1)),
  modal_clusters = vapply(chain_summaries, function(x) {
    as.integer(names(which.max(x$n_clusters)))
  }, integer(1))
)

## -----------------------------------------------------------------------------
cluster_count_trace <- function(fit) {
  vapply(fit$active_cluster_samples, sum, numeric(1))
}

constraint_trace <- function(fit) {
  vapply(fit$constraint_samples, constraint_to_model, character(1))
}

cluster_traces <- lapply(fits, cluster_count_trace)
constraint_traces <- lapply(fits, constraint_trace)

cluster_traces
constraint_traces

## ----fig.width = 7, fig.height = 4, fig.alt = "Trace plot of sampled cluster counts across two short chains."----
old_par <- par(mar = c(4, 4, 3, 1))
plot(
  cluster_traces[[1]],
  type = "b",
  pch = 19,
  ylim = range(unlist(cluster_traces)),
  col = "#0072B2",
  xlab = "Saved iteration",
  ylab = "Active clusters",
  main = "Cluster-count trace"
)
lines(cluster_traces[[2]], type = "b", pch = 19, col = "#D55E00")
legend("topright", legend = names(fits), col = c("#0072B2", "#D55E00"), lty = 1, pch = 19, bty = "n")
par(old_par)

## -----------------------------------------------------------------------------
co_clustering_matrix <- function(fit) {
  n <- length(fit$allocation_samples[[1]])
  out <- matrix(0, n, n)

  for (allocation in fit$allocation_samples) {
    out <- out + outer(allocation, allocation, "==")
  }

  out / length(fit$allocation_samples)
}

co_mat <- co_clustering_matrix(fits[[1]])
round(co_mat[1:6, 1:6], 2)

## ----fig.width = 5.5, fig.height = 5, fig.alt = "Heatmap of posterior co-clustering probabilities."----
image(
  seq_len(nrow(co_mat)),
  seq_len(ncol(co_mat)),
  co_mat[nrow(co_mat):1, ],
  col = hcl.colors(20, "YlGnBu", rev = TRUE),
  xlab = "Observation",
  ylab = "Observation",
  main = "Posterior co-clustering"
)

