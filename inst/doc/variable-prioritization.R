## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library(bpgmm)

simulate_screening_data <- function(n_per_cluster = 18, p = 6) {
  means <- rbind(
    c(-2.5, -1.5,  0.0, 0, 0, 0),
    c( 2.0, -0.5,  1.8, 0, 0, 0),
    c( 0.0,  2.2, -1.8, 0, 0, 0)
  )

  covariances <- list(
    diag(c(0.35, 0.35, 0.35, 1.40, 0.35, 1.20)),
    diag(c(0.35, 0.35, 0.35, 0.35, 1.40, 1.20)),
    matrix(c(
      0.35, 0,    0,    0,    0,    0,
      0,    0.35, 0,    0,    0,    0,
      0,    0,    0.35, 0,    0,    0,
      0,    0,    0,    1.00, 0.45, 0,
      0,    0,    0,    0.45, 1.00, 0,
      0,    0,    0,    0,    0,    1.20
    ), nrow = p)
  )

  n <- n_per_cluster * 3
  X <- matrix(NA_real_, nrow = p, ncol = n)
  true_cluster <- rep(seq_len(3), each = n_per_cluster)

  column <- 1
  for (k in seq_len(3)) {
    for (i in seq_len(n_per_cluster)) {
      X[, column] <- MASS::mvrnorm(1, mu = means[k, ], Sigma = covariances[[k]])
      column <- column + 1
    }
  }

  rownames(X) <- paste0("variable_", seq_len(p))
  list(X = X, true_cluster = true_cluster)
}

set.seed(2028)
sim <- simulate_screening_data()
X <- sim$X
true_cluster <- sim$true_cluster

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
    verbose = FALSE
  )
})
tail(fit_log, 1)

fit_summary <- summarize_pgmm_rjmcmc(fit, true_cluster = true_cluster)
fit_summary$ari

## -----------------------------------------------------------------------------
cluster_separation <- function(X, allocation) {
  vapply(seq_len(nrow(X)), function(j) {
    x <- X[j, ]
    overall <- mean(x)
    total <- sum((x - overall)^2)
    between <- sum(vapply(split(x, allocation), function(group) {
      length(group) * (mean(group) - overall)^2
    }, numeric(1)))
    if (total == 0) 0 else between / total
  }, numeric(1))
}

separation <- cluster_separation(X, fit_summary$allocation)
separation_table <- data.frame(
  variable = rownames(X),
  separation = round(separation, 3)
)
separation_table[order(separation_table$separation, decreasing = TRUE), ]

## ----fig.width = 6.5, fig.height = 4, fig.alt = "Bar plot of exploratory variable separation scores."----
ordered <- order(separation, decreasing = TRUE)
barplot(
  separation[ordered],
  names.arg = rownames(X)[ordered],
  las = 2,
  col = "#56B4E9",
  border = NA,
  ylab = "Between-cluster variation / total variation",
  main = "Exploratory variable separation"
)

## -----------------------------------------------------------------------------
loading_importance <- function(fit) {
  totals <- NULL
  count <- 0

  for (sample_id in seq_along(fit$lambda_samples)) {
    lambdas <- fit$lambda_samples[[sample_id]]
    active <- which(fit$active_cluster_samples[[sample_id]] == 1)

    for (k in active) {
      lambda <- lambdas[[k]]
      if (!is.matrix(lambda)) {
        next
      }
      score <- rowMeans(abs(lambda))
      if (is.null(totals)) {
        totals <- numeric(length(score))
      }
      totals <- totals + score
      count <- count + 1
    }
  }

  if (count == 0) {
    return(rep(NA_real_, nrow(X)))
  }
  totals / count
}

loading_score <- loading_importance(fit)
loading_table <- data.frame(
  variable = rownames(X),
  mean_abs_loading = round(loading_score, 3)
)
loading_table[order(loading_table$mean_abs_loading, decreasing = TRUE), ]

## ----fig.width = 6.5, fig.height = 4, fig.alt = "Bar plot of average posterior absolute loading magnitudes by variable."----
ordered_loading <- order(loading_score, decreasing = TRUE)
barplot(
  loading_score[ordered_loading],
  names.arg = rownames(X)[ordered_loading],
  las = 2,
  col = "#E69F00",
  border = NA,
  ylab = "Mean absolute loading",
  main = "Posterior loading magnitude"
)

