#' Hyperparameter set for the Bayesian PGMM sampler.
#'
#' @import methods
#' @slot alpha1 First Dirichlet prior parameter for component weights.
#' @slot alpha2 Second Dirichlet prior parameter for component weights.
#' @slot delta Shape parameter used in prior updates.
#' @slot ggamma Prior rate parameter used in covariance updates.
#' @slot bbeta Prior scale parameter used in covariance updates.
#' @keywords internal
#'
setClass(
  "Hparam",
  slots = c(
    alpha1 = "numeric",
    alpha2 = "numeric",
    delta  = "numeric",
    ggamma = "numeric",
    bbeta  = "numeric"
  ),
  prototype = list(
    alpha1 = numeric(),
    alpha2 = numeric(),
    delta  = numeric(),
    ggamma = numeric(),
    bbeta  = numeric()
  )
)


setValidity("Hparam", function(object) {
  if (object@alpha1 < 0 |
    object@alpha2 < 0 |
    object@delta < 0 |
    object@ggamma < 0 |
    object@bbeta < 0) {
    "Hyperparameter should be non-negative!"
  }
})
