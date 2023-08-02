#' Pad a Matrix
#'
#' @param x A matrix-like object
#' @param rowidx,colidx Integer indices for rows and columns from
#' \code{x} to \code{shape}
#' @param shape Target shape to pad out to
#' @template param-dots-method
#'
#' @return A matrix containing the data of \code{x} padded out to \code{shape}
#'
#' @noRd
pad_matrix <- function(x, ...) {
  UseMethod(generic = 'pad_matrix', object = x)
}

#' Pad a sparse Matrix with additional rows or columns
#'
#' @param x A dgTMatrix
#' @param colnames,rownames A vector of column or row names
#' to add to the matrix.
#' @param returns A padded matrix containing all provided
#' row/column names
#'
#' @importFrom Matrix sparseMatrix
#' @method pad_matrix default
#' @noRd
pad_matrix.default <- function(x, rownames = NULL, colnames = NULL, ...) {
  stopifnot(
    inherits(x, "Matrix"),
    is.character(colnames) || is.character(rownames)
  )
  # lookup table for Matrix representations
  mat_rep <- switch(
    EXPR = class(x),
    dgTMatrix = "T",
    dgCMatrix = "C",
    dgRMatrix = "R",
    stop("Untested Matrix object representation")

  )
  new_rownames <- setdiff(rownames, rownames(x))
  new_colnames <- setdiff(colnames, colnames(x))
  dtype <- typeof(methods::slot(object = x, name = 'x'))
  if (!is_empty(new_rownames)) {
    rpad <- Matrix::sparseMatrix(
      i = integer(0L),
      j = integer(0L),
      x = vector(mode = dtype, length = 0L),
      dims = c(length(new_rownames), ncol(x)),
      dimnames = list(new_rownames, colnames(x)),
      repr = mat_rep
    )
    x <- rbind(x, rpad)
  }
  if (!is_empty(new_colnames)) {
    cpad <- Matrix::sparseMatrix(
      i = integer(0L),
      j = integer(0L),
      x = vector(mode = dtype, length = 0L),
      dims = c(nrow(x), length(new_colnames)),
      dimnames = list(rownames(x), new_colnames),
      repr = mat_rep
    )
    x <- cbind(x, cpad)
  }
  x
}

#' @param sparse Return a \link[Matrix:TsparseMatrix-class]{sparse matrix}
#' @rdname pad_matrix
#' @method pad_matrix matrix
#' @noRd

pad_matrix.matrix <- function(x, rowidx, colidx, shape, sparse = FALSE, ...) {
  stopifnot(
    rlang::is_integerish(shape, n = 2L, finite = TRUE) && all(shape > 0L),
    rlang::is_integerish(rowidx, finite = TRUE) && all(rowidx > 0L),
    rlang::is_integerish(colidx, finite = TRUE) && all(colidx > 0L),
    all(dim(x) == c(length(rowidx), length(colidx))),
    is_scalar_logical(sparse)
  )
  if (!all(rowidx) <= shape[1L]) {
    stop('rowidx')
  } else if (!all(colidx <= shape[2L])) {
    stop('colidx')
  }
  type <- typeof(x)
  type <- match.arg(arg = type, choices = c('integer', 'double', 'logical'))
  mat <- if (isTRUE(sparse)) {
    Matrix::sparseMatrix(
      i = integer(),
      j = integer(),
      x = switch(EXPR = type, logical = logical(), numeric()),
      dims = shape,
      repr = 'T'
    )
  } else {
    matrix(
      data = vector(mode = type, length = prod(shape)),
      nrow = shape[1L],
      ncol = shape[2L]
    )
  }
  mat[rowidx, colidx] <- x
  return(mat)
}
