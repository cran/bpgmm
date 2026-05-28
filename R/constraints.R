#' Convert PGMM Constraint Codes to Paper Model Names
#'
#' The paper represents the eight PGMM covariance structures with three-letter
#' model names. Each letter is either `C` for constrained or `U` for
#' unconstrained. The letters indicate whether the loading matrix is shared
#' across clusters, whether the noise covariance is shared across clusters, and
#' whether the noise covariance is isotropic within clusters.
#'
#' @param constraint Integer or numeric vector of length three with entries
#'   `0` or `1`. `1` maps to `C`; `0` maps to `U`.
#'
#' @return A character scalar, one of `CCC`, `CCU`, `CUC`, `CUU`, `UCC`,
#'   `UCU`, `UUC`, or `UUU`.
#' @export
constraint_to_model <- function(constraint) {
  if (!is.numeric(constraint) && !is.integer(constraint)) {
    stop("constraint must be a numeric or integer vector", call. = FALSE)
  }
  if (length(constraint) != 3) {
    stop("constraint must have length 3", call. = FALSE)
  }
  if (any(is.na(constraint)) || any(!(constraint %in% c(0, 1)))) {
    stop("constraint entries must be 0/1 values", call. = FALSE)
  }

  paste0(ifelse(as.integer(constraint) == 1L, "C", "U"), collapse = "")
}

#' Convert PGMM Paper Model Names to Constraint Codes
#'
#' @param model Character scalar naming one of the eight PGMM covariance
#'   structures: `CCC`, `CCU`, `CUC`, `CUU`, `UCC`, `UCU`, `UUC`, or `UUU`.
#'
#' @return Integer vector of length three. `1` means constrained and `0` means
#'   unconstrained.
#' @export
model_to_constraint <- function(model) {
  valid_models <- c("CCC", "CCU", "CUC", "CUU", "UCC", "UCU", "UUC", "UUU")

  if (!is.character(model) || length(model) != 1L || is.na(model)) {
    stop("model must be a character scalar", call. = FALSE)
  }

  model <- toupper(model)
  if (!(model %in% valid_models)) {
    stop(
      "model must be one of: ",
      paste(valid_models, collapse = ", "),
      call. = FALSE
    )
  }

  as.integer(strsplit(chartr("CU", "10", model), "", fixed = TRUE)[[1]])
}
