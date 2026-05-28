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

## ----fig.width = 5.5, fig.height = 4, fig.alt = "Scatter plot of the small two-cluster example."----
plot(
  X[1, ], X[2, ],
  col = c("#0072B2", "#D55E00")[known_labels],
  pch = 19,
  xlab = "Variable 1",
  ylab = "Variable 2",
  main = "Small two-cluster example",
  asp = 1
)
legend(
  "topleft",
  legend = paste("Reference", sort(unique(known_labels))),
  col = c("#0072B2", "#D55E00"),
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
    constraint = model_to_constraint("UUU"),
    m_step = 0,
    v_step = 0,
    verbose = FALSE
  )
})
tail(fit_log, 1)

## -----------------------------------------------------------------------------
summary <- summarize_pgmm_rjmcmc(fit, true_cluster = known_labels)

summary$allocation
summary$n_clusters
summary$n_constraints
summary$ari

## ----fig.width = 5.5, fig.height = 4, fig.alt = "Scatter plot colored by posterior modal allocation."----
plot(
  X[1, ], X[2, ],
  col = c("#009E73", "#CC79A7", "#E69F00")[summary$allocation],
  pch = 19,
  xlab = "Variable 1",
  ylab = "Variable 2",
  main = "Posterior modal allocation",
  asp = 1
)
text(X[1, ], X[2, ], labels = seq_along(summary$allocation), pos = 3, cex = 0.75)
legend(
  "topleft",
  legend = paste("Cluster", sort(unique(summary$allocation))),
  col = c("#009E73", "#CC79A7", "#E69F00")[sort(unique(summary$allocation))],
  pch = 19,
  bty = "n"
)

## -----------------------------------------------------------------------------
citation("bpgmm")

