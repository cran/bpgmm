## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library(bpgmm)

set.seed(2026)

X <- cbind(
  matrix(rnorm(8, mean = -2, sd = 0.2), nrow = 2),
  matrix(rnorm(8, mean = 2, sd = 0.2), nrow = 2)
)
known_labels <- rep(1:2, each = 4)

dim(X)
known_labels

## ----fig.width = 5.5, fig.height = 4, fig.alt = "Scatter plot of simulated observations colored by known cluster."----
cluster_cols <- c("#0072B2", "#D55E00", "#009E73")
plot(
  X[1, ], X[2, ],
  col = cluster_cols[known_labels],
  pch = 19,
  xlab = "Variable 1",
  ylab = "Variable 2",
  main = "Simulated data",
  asp = 1
)
legend(
  "topleft",
  legend = paste("Known", sort(unique(known_labels))),
  col = cluster_cols[sort(unique(known_labels))],
  pch = 19,
  bty = "n"
)

## -----------------------------------------------------------------------------
fit_log <- capture.output({
  fit <- pgmm_rjmcmc(
    X = X,
    m_init = 2,
    m_range = c(1, 3),
    q_new = 1,
    burn = 1,
    niter = 3,
    constraint = "UUU",
    m_step = 0,
    v_step = 0,
    verbose = FALSE
  )
})
tail(fit_log, 1)

## -----------------------------------------------------------------------------
names(fit)
length(fit$allocation_samples)

## -----------------------------------------------------------------------------
fit_summary <- summarize_pgmm_rjmcmc(fit, true_cluster = known_labels)

fit_summary$allocation
fit_summary$n_clusters
fit_summary$n_constraints
fit_summary$ari

## ----fig.width = 7, fig.height = 4, fig.alt = "Two-panel plot comparing known labels and posterior allocation."----
old_par <- par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
plot(
  X[1, ], X[2, ],
  col = cluster_cols[known_labels],
  pch = 19,
  xlab = "Variable 1",
  ylab = "Variable 2",
  main = "Known labels",
  asp = 1
)
plot(
  X[1, ], X[2, ],
  col = cluster_cols[fit_summary$allocation],
  pch = 19,
  xlab = "Variable 1",
  ylab = "Variable 2",
  main = "Posterior allocation",
  asp = 1
)
par(old_par)

## -----------------------------------------------------------------------------
last <- length(fit$tau_samples)

fit$tau_samples[[last]]
fit$mean_samples[[last]]
lapply(fit$lambda_samples[[last]], dim)
lapply(fit$psi_samples[[last]], dim)

## ----eval = FALSE-------------------------------------------------------------
# fit <- pgmm_rjmcmc(
#   X = your_data_matrix,
#   m_init = 2,
#   m_range = c(1, 10),
#   q_new = 4,
#   burn = 5000,
#   niter = 15000,
#   constraint = model_to_constraint("UUU"),
#   m_step = 1,
#   v_step = 1,
#   split_combine = 1,
#   verbose = FALSE
# )
# 
# summarize_pgmm_rjmcmc(fit, true_cluster = known_labels)

