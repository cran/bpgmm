## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library(bpgmm)

iris_numeric <- as.matrix(iris[, 1:4])
iris_labels <- as.integer(iris$Species)

dim(iris_numeric)

X <- t(iris_numeric)
dim(X)

## -----------------------------------------------------------------------------
rownames(X)
X[, 1:3]

## -----------------------------------------------------------------------------
all(is.finite(X))
storage.mode(X)

## -----------------------------------------------------------------------------
X_scaled <- t(scale(t(X)))

round(rowMeans(X_scaled), 6)
round(apply(X_scaled, 1, sd), 6)

## -----------------------------------------------------------------------------
p <- nrow(X_scaled)
q_new <- min(2, p - 1)
q_new

## -----------------------------------------------------------------------------
m_init <- 3
m_range <- c(1, 5)

## ----fig.width = 5.5, fig.height = 4.2, fig.alt = "Scatter plot of scaled iris variables colored by species."----
species_cols <- c("#0072B2", "#D55E00", "#009E73")
plot(
  X_scaled[1, ], X_scaled[2, ],
  col = species_cols[iris_labels],
  pch = 19,
  xlab = rownames(X_scaled)[1],
  ylab = rownames(X_scaled)[2],
  main = "Scaled iris data",
  asp = 1
)
legend(
  "topleft",
  legend = levels(iris$Species),
  col = species_cols,
  pch = 19,
  bty = "n"
)

## -----------------------------------------------------------------------------
model_to_constraint("UUU")
constraint_to_model(c(1, 1, 1))

## ----eval = FALSE-------------------------------------------------------------
# fit <- pgmm_rjmcmc(
#   X = X_scaled,
#   m_init = m_init,
#   m_range = m_range,
#   q_new = q_new,
#   burn = 1000,
#   niter = 5000,
#   constraint = "UUU",
#   m_step = 1,
#   v_step = 1,
#   split_combine = 1,
#   verbose = FALSE
# )
# 
# summarize_pgmm_rjmcmc(fit, true_cluster = iris_labels)

