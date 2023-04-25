.anndata_to_sce_rd <- function(x) {
  if (is.null(x)) {
    return(NULL)
  }
  stopifnot(is.character(x))
  return(toupper(gsub(pattern = '^X_', replacement = '', x = x)))
}

.check_sce_installed <- function(quietly = FALSE) {
  pkg <- 'SingleCellExperiment'
  checks <- c(
    installed = requireNamespace(pkg, quietly = TRUE),
    version = tryCatch(
      expr = utils::packageVersion(pkg) >= .MINIMUM_SCE_VERSION(),
      error = function(...) {
        return(FALSE)
      }
    )
  )
  if (isTRUE(quietly)) {
    return(invisible(all(checks)))
  }
  if (!checks['installed']) {
    stop(sQuote(pkg), " must be installed", call. = FALSE)
  }
  if (!checks['version']) {
    stop(
      sQuote(pkg),
      " must be version ",
      .MINIMUM_SCE_VERSION('c'),
      " or higher",
      call. = FALSE
    )
  }
  return(invisible(TRUE))
}

.MINIMUM_SCE_VERSION <- function(repr = c('v', 'c')) {
  repr <- repr[1L]
  repr <- match.arg(repr)
  version <- '1.20.0'
  return(switch(EXPR = repr, v = package_version(version), version))
}
