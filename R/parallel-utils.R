available_cores <- function() {
  cores <- parallel::detectCores(logical = FALSE)
  if (is.na(cores) || cores < 1L) {
    return(1L)
  }
  as.integer(cores)
}

validate_positive_integer <- function(x, name) {
  if (!isTRUE(length(x) == 1L) || !is.finite(x) || x < 1) {
    stop(name, " must be a positive scalar", call. = FALSE)
  }
  as.integer(x)
}

restore_random_seed <- function(seed) {
  if (is.null(seed)) {
    if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
    return(invisible(NULL))
  }
  assign(".Random.seed", seed, envir = .GlobalEnv)
  invisible(NULL)
}
