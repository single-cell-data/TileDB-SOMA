#' Convert a sparse matrix to a COO-formatted Data Frame
#'
#' A coordinate-oriented (COO) data frame is similar to the *triplet
#' representation* of a sparse matrix, but string dimensions (if present) are
#' maintained in the `i`/`j` columns. In addition, a list of multiple sparse
#' matrices can be consolidated into a single COO data frame if they are
#' *layerable* (i.e., if they have the same dimensions and non-zero
#' coordinates).
#'
#' @param x any matrix-like object coercible to a `TsparseMatrix` or a list of
#' matrix-like objects that all share the same dimensions and non-zero
#' coordinates (i.e., are *layerable*).
#' @returns A `data.frame` with factor columns for the i/j indices, and a value
#' column for each of the matrices included in `x`
#' @importFrom Matrix mat2triplet
#' @noRd
matrix_to_coo <- function(x, index_cols = c("i", "j"), value_cols = NULL) {
  if (!is.list(x)) x <- list(x)

  if (is.null(value_cols)) {
    if (is.null(names(x))) {
      value_cols <- paste0("value", seq_along(x))
    } else {
      value_cols <- names(x)
    }
  }

  stopifnot(length(index_cols) == 2, length(value_cols) == length(x))
  is_a_matrix <- vapply_lgl(x, is_matrix)
  if (!all(is_a_matrix)) {
    stop("When 'x' is a list all elements must contain a matrix-like object")
  }

  # Coerce all matrices to TsparseMatrix so they share a common coordinate space
  # when testing for layerability. The same coercion is performed within
  # mat2triplet so we're just paying the computational cost here instead.
  x <- lapply(x, as, Class = "TsparseMatrix")

  # Convert all matrices to COO format.
  # Add 1 column of values to the coo df for each additional layerable matrix
  for (i in seq_along(x)) {
    value_col <- value_cols[i]
    if (!are_layerable(x[[1]], x[[i]])) {
      stop(sprintf("Matrix 1 and %i are not layerable", i))
    }

    # convert to COO format w/ user specified index and value column names
    coo_i <- as.data.frame(
      setNames(
        Matrix::mat2triplet(x[[i]]),
        nm = c(index_cols, value_col)
      )
    )


    # When there are multiple matrices we need to order coordinates by index
    # values so we can efficiently/correctly merge them together
    if (length(x) > 1) {
      coord_order <- order(
        coo_i[[index_cols[1]]],
        coo_i[[index_cols[2]]]
      )
      coo_i <- coo_i[coord_order,]
    }

    if (i == 1) {
      coo <- coo_i
    } else {
      coo[[value_col]] <- coo_i[[value_col]]
    }
  }

  # Dimension names from the first matrix are used for the coo index columns.

  # To avoid materializing strings values for the index columns we use factors
  # instead. We don't need all the overhead of factor() because we already have
  # integer indices that map to the dimension names of the matrix, so we build
  # them manually using structure().

  i_col <- index_cols[1]
  j_col <- index_cols[2]

  coo[[i_col]] <- structure(
    coo[[i_col]],
    levels = rownames(x[[1]]),
    class = "factor"
  )
  coo[[j_col]] <- structure(
    coo[[j_col]],
    levels = colnames(x[[1]]),
    class = "factor"
  )

  # reset rownames in case they were ordered
  as.data.frame(coo, row.names = NULL)
}
