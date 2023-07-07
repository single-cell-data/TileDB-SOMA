.anndata_to_sce_reduced_dim <- function(x) {
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

#' Matrix to \code{S4Vectors::SelfHits}
#'
#' Re-implement \pkg{SingleCellExperiment}'s internal \code{.mat2hits}
#' function as theirs doesn't recognize the difference between \code{Matrix}
#' and \code{matrix} without \pkg{Matrix} being attached. Also correct cases
#' where \code{S4Vectors::SelfHits(from =} and \code{S4Vectors::SelfHits(to =}
#' are not integers by explicitly casting to integer
#'
#' @param mat A matrix
#'
#' @return An \code{\link[S4Vectors:SelfHits]{S4Vectors::SelfHits}} object
#'
#' @keywords internal
#'
#' @noRd
#'
.mat_to_hits <- function(mat) {
  stopifnot(
    "S4Vectors must be installed" = requireNamespace('S4Vectors', quietly = TRUE),
    "'mat' must be a matrix" = is_matrix(mat)
  )
  f <- if (inherits(mat, 'Matrix')) {
    Matrix::which
  } else {
    base::which
  }
  i <- f(mat != 0, arr.ind = TRUE)
  return(S4Vectors::SelfHits(
    from = as.integer(i[, 1L]),
    to = as.integer(i[, 2L]),
    nnode = nrow(mat),
    x = mat[i]
  ))
}

.MINIMUM_SCE_VERSION <- function(repr = c('v', 'c')) {
  repr <- repr[1L]
  repr <- match.arg(repr)
  version <- '1.20.0'
  return(switch(EXPR = repr, v = package_version(version), version))
}
