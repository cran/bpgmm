#' (internal)
#' @noRd
.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    paste0(
      "bpgmm ",
      utils::packageVersion("bpgmm"),
      " loaded. If you use bpgmm in published work, please cite it with citation(\"bpgmm\")."
    )
  )
}
