# matrixZeroBasedView is a wrapper shim for a matrix or Matrix::sparseMatrix providing
# elemental access using zero-based indexes.
#
# For M0 <- matrixZeroBasedView(M1):
# - M0[i,j] is equivalent to M1[i+1,j+1].
# - In particular, if M0[i,j] is vector- or matrix-valued, then the returned
#   vector/matrix is ONE-based.
# - as.one.based(M0) returns M1.
# - The only other supported operations are: dim(M0), nrow(M0), ncol(M0).

#' Zero-based matrix shim
#'
#' @param one_based_matrix matrix or Matrix::sparseMatrix
#'
#' @return Shim providing elemental access to the matrix using zero-based indexes.
#' @export
#'
#' @examples
matrixZeroBasedView <- function(one_based_matrix) {
  if (!inherits(one_based_matrix, "matrix") && !inherits(one_based_matrix, "sparseMatrix")) {
    stop("Matrix object must inherit class matrix or Matrix::sparseMatrix.")
  }
  structure(list(one_based_matrix = one_based_matrix), class = "matrixZeroBasedView")
}

#' Zero-based matrix element access
#'
#' @param x
#' @param i
#' @param j
#' @param drop
#'
#' @return
#' @export
#'
#' @examples
`[.matrixZeroBasedView` <- function(x, i, j, drop = TRUE) {
  if (missing(i) && missing(j)) {
    x$one_based_matrix[, , drop = drop]
  } else if (missing(i)) {
    x$one_based_matrix[, j + 1, drop = drop]
  } else if (missing(j)) {
    x$one_based_matrix[i + 1, , drop = drop]
  } else {
    x$one_based_matrix[i + 1, j + 1, drop = drop]
  }
}

#' dim
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
dim.matrixZeroBasedView <- function(x) {
  dim(x$one_based_matrix)
}

#' nrow
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
nrow.matrixZeroBasedView <- function(x) {
  nrow(x$one_based_matrix)
}

#' ncol
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
ncol.matrixZeroBasedView <- function(x) {
  ncol(x$one_based_matrix)
}

#' Get one-based object
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
as.one.based <- function(x) {
  UseMethod("as.one.based")
}

#' Get one-based matrix object underlying zero-based view
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
as.one.based.matrixZeroBasedView <- function(x) {
  x$one_based_matrix
}
