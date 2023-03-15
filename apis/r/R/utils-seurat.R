.check_seurat_installed <- function(quietly = FALSE) {
  pkg <- 'SeuratObject'
  checks <- c(
    installed = requireNamespace(pkg, quietly = TRUE),
    version = tryCatch(
      expr = utils::packageVersion(pkg) >= .MINIMUM_SEURAT_VERSION(),
      error = \(...) FALSE
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
      .MINIMUM_SEURAT_VERSION('c'),
      " or higher",
      call. = FALSE
    )
  }
  return(invisible(TRUE))
}

.MINIMUM_SEURAT_VERSION <- function(repr = c('v', 'c')) {
  repr <- repr[1L]
  repr <- match.arg(arg = repr)
  version <- '4.1.0'
  return(switch(EXPR = repr, v = package_version(version), version))
}
