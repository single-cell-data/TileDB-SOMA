.anndata_to_sce_reduced_dim <- function(x) {
  if (is.null(x)) {
    return(NULL)
  }
  stopifnot(is.character(x))
  return(toupper(gsub(pattern = '^X_', replacement = '', x = x)))
}

#' \code{S4Vectors::Hits} to Matrix
#'
#' Re-implement \pkg{SingleCellExperiment}'s internal \code{.hits2mat}
#' function. This version returns a \code{TsparseMatrix} instead of
#' a \code{CsparseMatrix}
#'
#' @param hits A \code{Hits} object
#'
#' @return A \code{TsparseMatrix} representation of \code{hits}
#'
#' @keywords internal
#'
#' @noRd
#'
.hits_to_mat <- function(hits) {
  stopifnot(
    "S4Vectors must be installed" = requireNamespace('S4Vectors', quietly = TRUE),
    "'hits' must be a 'Hits' object" = inherits(hits, 'Hits')
  )
  meta_cols <- S4Vectors::mcols(hits)
  x <- if (ncol(meta_cols)) {
    meta_cols[[1L]]
  } else {
    rep.int(TRUE, times = length(hits))
  }
  return(Matrix::sparseMatrix(
    i = S4Vectors::queryHits(hits),
    j = S4Vectors::subjectHits(hits),
    x = x,
    repr = 'T',
    dims = rep.int(S4Vectors::nnode(hits), times = 2L),
    use.last.ij = TRUE
  ))
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
