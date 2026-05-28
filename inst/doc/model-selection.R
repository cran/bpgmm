## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library(bpgmm)

simulate_mfa_data <- function(n_per_cluster = 20, p = 6, q = 2) {
  means <- rbind(
    c(-3.0, -2.0,  0.0, 0, 0, 0),
    c( 2.5, -1.0,  2.0, 0, 0, 0),
    c( 0.0,  2.5, -2.0, 0, 0, 0)
  )

  lambdas <- list(
    matrix(c(1.2,  0.5, 0.2, 0, 0, 0,
             0.0,  0.2, 0.8, 0.3, 0, 0), nrow = p),
    matrix(c(0.2,  1.1, 0.6, 0, 0, 0,
             0.8,  0.1, 0.2, 0.6, 0, 0), nrow = p),
    matrix(c(0.7, -0.4, 1.0, 0, 0, 0,
            -0.2,  0.8, 0.4, 0.3, 0, 0), nrow = p)
  )

  psi <- list(
    diag(c(0.25, 0.35, 0.30, 0.80, 0.90, 1.00)),
    diag(c(0.30, 0.25, 0.40, 0.80, 0.90, 1.00)),
    diag(c(0.35, 0.30, 0.25, 0.80, 0.90, 1.00))
  )

  n <- n_per_cluster * 3
  X <- matrix(NA_real_, nrow = p, ncol = n)
  true_cluster <- rep(seq_len(3), each = n_per_cluster)

  column <- 1
  for (k in seq_len(3)) {
    for (i in seq_len(n_per_cluster)) {
      latent_score <- rnorm(q)
      noise <- MASS::mvrnorm(1, mu = rep(0, p), Sigma = psi[[k]])
      X[, column] <- means[k, ] + lambdas[[k]] %*% latent_score + noise
      column <- column + 1
    }
  }

  rownames(X) <- paste0("variable_", seq_len(p))
  list(X = X, true_cluster = true_cluster)
}

set.seed(2027)
sim <- simulate_mfa_data()
X <- sim$X
true_cluster <- sim$true_cluster

dim(X)
table(true_cluster)

## ----fig.width = 5.5, fig.height = 4.2, fig.alt = "Scatter plot of the first two variables colored by true cluster."----
cluster_cols <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00")
plot(
  X[1, ], X[2, ],
  col = cluster_cols[true_cluster],
  pch = 19,
  xlab = rownames(X)[1],
  ylab = rownames(X)[2],
  main = "Simulated MFA data",
  asp = 1
)
legend(
  "topleft",
  legend = paste("True cluster", sort(unique(true_cluster))),
  col = cluster_cols[sort(unique(true_cluster))],
  pch = 19,
  bty = "n"
)

## -----------------------------------------------------------------------------
fit_log <- capture.output({
  fit <- pgmm_rjmcmc(
    X = X,
    m_init = 3,
    m_range = c(1, 5),
    q_new = 2,
    burn = 2,
    niter = 6,
    constraint = "UUU",
    m_step = 1,
    v_step = 1,
    split_combine = 0,
    verbose = FALSE
  )
})
tail(fit_log, 1)

## -----------------------------------------------------------------------------
fit_summary <- summarize_pgmm_rjmcmc(fit, true_cluster = true_cluster)

fit_summary$n_clusters
fit_summary$n_constraints
fit_summary$ari

## ----fig.width = 7, fig.height = 4, fig.alt = "Bar plots of posterior counts for cluster number and covariance model."----
old_par <- par(mfrow = c(1, 2), mar = c(5, 4, 3, 1))
barplot(
  fit_summary$n_clusters,
  col = "#56B4E9",
  border = NA,
  xlab = "Number of clusters",
  ylab = "Posterior sample count",
  main = "Cluster-number selection"
)
barplot(
  fit_summary$n_constraints,
  col = "#E69F00",
  border = NA,
  xlab = "PGMM model",
  ylab = "Posterior sample count",
  main = "Covariance-model selection"
)
par(old_par)

## ----fig.width = 7, fig.height = 4, fig.alt = "Two-panel scatter plot comparing true cluster labels and posterior modal allocation."----
old_par <- par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
plot(
  X[1, ], X[2, ],
  col = cluster_cols[true_cluster],
  pch = 19,
  xlab = rownames(X)[1],
  ylab = rownames(X)[2],
  main = "True clusters",
  asp = 1
)
plot(
  X[1, ], X[2, ],
  col = cluster_cols[fit_summary$allocation],
  pch = 19,
  xlab = rownames(X)[1],
  ylab = rownames(X)[2],
  main = "Posterior modal allocation",
  asp = 1
)
par(old_par)

