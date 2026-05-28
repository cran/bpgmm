## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----fig.width = 5.5, fig.height = 4.5, fig.alt = "Two covariance ellipses illustrating cluster-specific covariance matrices."----
ellipse_points <- function(center, sigma, level = 0.90, n = 120) {
  eig <- eigen(sigma, symmetric = TRUE)
  angles <- seq(0, 2 * pi, length.out = n)
  circle <- rbind(cos(angles), sin(angles))
  radius <- sqrt(qchisq(level, df = 2))
  points <- t(center + radius * eig$vectors %*% diag(sqrt(eig$values), 2) %*% circle)
  colnames(points) <- c("x", "y")
  points
}

lambda_1 <- matrix(c(1.1, 0.4), ncol = 1)
lambda_2 <- matrix(c(0.2, 1.0), ncol = 1)
psi_1 <- diag(c(0.20, 0.08))
psi_2 <- diag(c(0.10, 0.25))
sigma_1 <- lambda_1 %*% t(lambda_1) + psi_1
sigma_2 <- lambda_2 %*% t(lambda_2) + psi_2

ell_1 <- ellipse_points(c(-1.2, -0.5), sigma_1)
ell_2 <- ellipse_points(c(1.1, 0.6), sigma_2)

plot(
  ell_1,
  type = "l",
  lwd = 2,
  col = "#0072B2",
  xlim = c(-3.5, 3.5),
  ylim = c(-2.5, 3),
  xlab = "Variable 1",
  ylab = "Variable 2",
  main = "Mixture-of-factor-analyzers covariance"
)
lines(ell_2, lwd = 2, col = "#D55E00")
points(rbind(c(-1.2, -0.5), c(1.1, 0.6)), pch = 19, col = c("#0072B2", "#D55E00"))
arrows(-1.2, -0.5, -1.2 + lambda_1[1], -0.5 + lambda_1[2], col = "#0072B2", lwd = 2, length = 0.08)
arrows(1.1, 0.6, 1.1 + lambda_2[1], 0.6 + lambda_2[2], col = "#D55E00", lwd = 2, length = 0.08)
legend(
  "topleft",
  legend = c("Cluster 1 covariance", "Cluster 2 covariance", "Loading direction"),
  col = c("#0072B2", "#D55E00", "gray30"),
  lwd = c(2, 2, 2),
  bty = "n"
)

## -----------------------------------------------------------------------------
library(bpgmm)

models <- c("CCC", "CCU", "CUC", "CUU", "UCC", "UCU", "UUC", "UUU")
data.frame(
  model = models,
  constraint = vapply(
    models,
    function(x) paste(model_to_constraint(x), collapse = ","),
    character(1)
  )
)

## -----------------------------------------------------------------------------
model_to_constraint("UUU")
constraint_to_model(c(1, 1, 0))

## ----eval = FALSE-------------------------------------------------------------
# fit <- pgmm_rjmcmc(
#   X = your_data_matrix,
#   m_init = 5,
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

## -----------------------------------------------------------------------------
citation("bpgmm")

