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
matrixZeroBasedView <- function(one_based_matrix) {
  if (!inherits(one_based_matrix, "matrix") && !inherits(one_based_matrix, "sparseMatrix")) {
    stop("Matrix object must inherit class matrix or Matrix::sparseMatrix.")
  }
  structure(list(one_based_matrix = one_based_matrix), class = "matrixZeroBasedView")
}

#' Zero-based matrix element access
#'
#' @param x The zero-based matrix view.
#' @param i Row index (zero-based).
#' @param j Column index (zero-based).
#' @param drop Coerce result to lowest possible dimension.
#'
#' @return The specified matrix elements. Vector- or matrix-valued results are returned
#'         as conventional one-based R objects.
#' @export
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
#' @param x The zero-based matrix view.
#'
#' @return The dimensions of the matrix.
#' @export
dim.matrixZeroBasedView <- function(x) {
  dim(x$one_based_matrix)
}

#' nrow
#'
#' @param x The zero-based matrix view.
#'
#' @return Matrix row count.
#' @export
nrow.matrixZeroBasedView <- function(x) {
  nrow(x$one_based_matrix)
}

#' ncol
#'
#' @param x The zero-based matrix view.
#'
#' @return Matrix column count.
#' @export
ncol.matrixZeroBasedView <- function(x) {
  ncol(x$one_based_matrix)
}

#' Get one-based object
#'
#' @param x The object.
#'
#' @return A one-based version/view of the object.
#' @export
as.one.based <- function(x) {
  UseMethod("as.one.based")
}

#' Get one-based matrix object underlying zero-based view
#'
#' @param x The zero-based matrix view.
#'
#' @return The one-based `matrix` or `Matrix::sparseMatrix` underlying the view.
#' @export
as.one.based.matrixZeroBasedView <- function(x) {
  x$one_based_matrix
}

#' Zero-based matrix element assigment
#'
#' This method errors as a read-only view is implemented.
#'
#' @param x The zero-based matrix view.
#' @param i Row index (zero-based).
#' @param j Column index (zero-based).
#' @param val The to-be assigned value.
#'
#' @return Nothing as the method errors.
#' @export
`[<-.matrixZeroBasedView` <- function(x, i, j, value) {
  stop("matrixZeroBasedView is read-only; use as.one.based() to get full-featured matrix object", call. = FALSE)
}
