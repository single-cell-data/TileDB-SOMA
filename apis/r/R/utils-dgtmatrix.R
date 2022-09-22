#' Convert from COO-formatted Data Frame to dgTMatrix
#' @param x A COO-formatted `data.frame` with columns for the i/j indices, and
#' and one or more value columns.
#' @returns A list of `dgTMatrix` objects, with one element per value column in
#' `x`.
#' @noRd
dataframe_to_dgtmatrix <- function(x, index_cols = c("i", "j")) {
  stopifnot(
    is.data.frame(x),
    length(index_cols) == 2,
    all(index_cols %in% colnames(x))
  )

  value_cols <- setdiff(colnames(x), index_cols)
  dim_labels <- as.list(x[index_cols])
  dim_names <- lapply(dim_labels, unique)
  dim_lengths <- vapply(dim_names, length, FUN.VALUE = integer(1L))

  mapply(
    FUN = Matrix::sparseMatrix,
    x = x[value_cols],
    MoreArgs = list(
      i = match(dim_labels[[1]], dim_names[[1]]),
      j = match(dim_labels[[2]], dim_names[[2]]),
      dims = dim_lengths,
      dimnames = unname(dim_names),
      repr = "T"
    )
  )
}

# Matrices with identical dimension names and non-empty coordinates can be
# stored as different layers (i.e., attributes of the same array)
#' @importFrom Matrix nnzero
are_layerable <- function(x, y) {
  stopifnot(is_matrix(x) && is_matrix(y))
  if (identical(x, y)) return(TRUE)
  dimnames_match <- identical(dimnames(x), dimnames(y))
  nonemptycells_match <- Matrix::nnzero(x) == Matrix::nnzero(y)
  dimnames_match && nonemptycells_match
}


#' Pad a sparse Matrix with additional rows or columns
#' @param x A dgTMatrix
#' @param colnames,rownames A vector of column or row names to add to the
#' matrix.
#' @param returns A padded matrix containing all provided row/column names
#' @importFrom Matrix sparseMatrix
#' @noRd
pad_matrix <- function(x, rownames = NULL, colnames = NULL) {
  stopifnot(
    inherits(x, "Matrix"),
    is.character(colnames) || is.character(rownames)
  )

  # lookup table for Matrix representations
  mat_rep <- switch(class(x),
    dgTMatrix = "T",
    dgCMatrix = "C",
    dgRMatrix = "R",
    stop("Untested Matrix object representation")
  )

  new_rownames <- setdiff(rownames, rownames(x))
  new_colnames <- setdiff(colnames, colnames(x))
  dtype <- typeof(x@x)

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
