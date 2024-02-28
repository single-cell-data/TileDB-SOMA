#' Transformer function: SOMAArray to Arrow table
#'
#' @description Converts the results of a \link{soma_array_reader} or
#' \link{sr_next} to an arrow::\link[arrow]{Table}
#' @param x A List object with two pointers to Arrow array data and schema
#' @return arrow::\link[arrow]{Table}
#' @noRd
soma_array_to_arrow_table <- function(x) {
  check_arrow_pointers(x)
  arrow::as_arrow_table(
    arrow::RecordBatch$import_from_c(x$array_data, x$schema)
  )
}

soma_array_to_arrow_table_concat <- function(it) {
  stopifnot(inherits(it, 'ReadIter'))
  tbl <- it$read_next()
  while (!it$read_complete()) {
    nxt <- it$read_next()
    tbl <- arrow::concat_tables(tbl, nxt)
  }
  return(tbl)
}

soma_array_to_sparse_matrix_concat <- function(obj, zero_based = FALSE) {
  if (obj$read_complete()) {
    warning("Iteration complete, returning NULL")
    return(NULL)
  }

  mat <- obj$read_next()

  while (!obj$read_complete()) {
    if (zero_based) {
      mat <- mat$sum(obj$read_next())
    } else {
      mat <- mat + obj$read_next()
    }
  }

  mat
}


#' Transformer function: Arrow table to Matrix::sparseMatrix
#'
#' @description Converts a \link[arrow]{Table} of sparse format (columns:
#' "soma_dim_0", "soma_dim_1", "soma_data") to a \link{matrixZeroBasedView}
#' @param tbl \link[arrow]{Table} with columns "soma_dim_0", "soma_dim_1", and
#' "soma_data"
#' @param repr Optional one-character code for sparse matrix representation type
#' @param shape Numerical vector with two elements, one for each dimension. If
#' \code{NULL}, then the following is used \code{1 + c(max(tbl["soma_dim_0"]),
#' max(tbl["soma_dim_1"]))}
#' @param repr Optional one-character code for sparse matrix representation type
#' @param zero_based Logical, if TRUE returns a Matrix::\link{sparseMatrix}
#' otherwise \link{matrixZeroBasedView} of Matrix::\link[Matrix]{sparseMatrix}
#' @return Matrix::\link{sparseMatrix} or \link{matrixZeroBasedView} of
#' Matrix::\link[Matrix]{sparseMatrix}
#' @noRd
arrow_table_to_sparse <- function(tbl, repr = c("C", "T", "R"), shape = NULL, zero_based = FALSE) {

  # To instantiate the one-based Matrix::sparseMatrix, we need to add 1 to the
  # zero-based soma_dim_0 and soma_dim_1 (done by arrow_table_to_sparse). But, because these dimensions are
  # usually populated with soma_joinid, users will need to access the matrix
  # using the original, possibly-zero IDs. Therefore, we'll wrap the one-based
  # sparseMatrix with a shim providing basic access with zero-based indexes.
  # If needed, user can then explicitly ask the shim for the underlying
  # sparseMatrix using `as.one.based()`.

  if (is.null(shape)) {
    shape <- c(max(tbl$soma_dim_0)$as_vector(), max(tbl$soma_dim_1)$as_vector())
  }

  stopifnot(
    "'shape' must not exceed '.Machine$integer.max'." =
      all(shape <= .Machine$integer.max),
    "A Matrix::sparseMatrix cannot hold more than 2^31 - 1 non-zero values" =
      nrow(tbl) <= .Machine$integer.max
  )

  exceedsInt32Limit <- (
    tbl$soma_dim_0 >= .Machine$integer.max | tbl$soma_dim_1 >= .Machine$integer.max
  )

  if (any(exceedsInt32Limit)$as_vector()) {
    stop(
      "Query contains 0-based coordinates outside '[0, 2^31 - 1)'.\n",
      "  - Matrix::sparseMatrix cannot be created with these values",
      call. = FALSE
    )
  }

  mat <- Matrix::sparseMatrix(i = tbl$soma_dim_0$as_vector(),
                              j = tbl$soma_dim_1$as_vector(),
                              x = tbl$soma_data$as_vector(),
                              dims = shape,
                              repr = repr,
                              index1 = FALSE)
  if (zero_based) {
      matrixZeroBasedView$new(mat)
  } else {
      mat
  }
}


#' Transformer function: Arrow table to matrix
#'
#' @description Converts a \link[arrow]{Table} of sparse format (columns:
#' "soma_dim_0", "soma_dim_1", "soma_data") to a \link{matrixZeroBasedView} of a
#' \link{matrix}.
#' @param tbl  \link[arrow]{Table} with columns "soma_dim_0", "soma_dim_1", and
#' "soma_data"
#' @param byrow Logical, TRUE if "soma_data" is ordered by row, this argument is
#' directly passed to the argument \code{byrow} of \link{matrix}
#' @return \link{matrixZeroBasedView} of \link[base]{matrix}
#' @noRd
arrow_table_to_dense <- function(tbl, byrow) {

  nrows <- length(unique(as.numeric(tbl$GetColumnByName("soma_dim_0"))))
  ncols <- length(unique(as.numeric(tbl$GetColumnByName("soma_dim_1"))))
  soma_data <- as.numeric(tbl$GetColumnByName("soma_data"))

  mat <- matrix(soma_data, nrow = nrows, ncol = ncols, byrow = byrow)
  matrixZeroBasedView$new(mat)
}
