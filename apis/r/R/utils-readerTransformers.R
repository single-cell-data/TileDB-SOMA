#' @description Converts the results of a \link{soma_array_reader} or \link{sr_next} to 
#' an \link[Arrow]{Table}
#' @return \link[Arrow]{Table}
soma_array_to_arrow_table <- function(x) {
    arrow::as_arrow_table(arrow::RecordBatch$import_from_c(x[[1]], x[[2]]))
}

#' @description Converts a \link[Arrow]{Table} of sparse format (columns: "soma_dim_0", 
#' "soma_dim_1", "soma_data") to a \link{matrixZeroBasedView}
#' @param tbl  \link[Arrow]{Table} with columns "soma_dim_0", "soma_dim_1", and "soma_datah"
#' @param repr Optional one-character code for sparse matrix representation type
#' @param dims_one_based Numerical vectors with two elements, one for each dimension. If
#' \code{NULL}, then the following is used \code{c(max(tbl["soma_dim_0"]), max(tbl["soma_dim_1"]))}
#' @return \link{matrixZeroBasedView}
arrow_table_to_sparse <- function(tbl, repr = c("C", "T", "R"), dims_one_based = NULL) {

  # To instantiate the one-based Matrix::sparseMatrix, we need to add 1 to the
  # zero-based soma_dim_0 and soma_dim_1 (done by arrow_table_to_sparse). But, because these dimensions are
  # usually populated with soma_joinid, users will need to access the matrix
  # using the original, possibly-zero IDs. Therefore, we'll wrap the one-based
  # sparseMatrix with a shim providing basic access with zero-based indexes.
  # If needed, user can then explicitly ask the shim for the underlying
  # sparseMatrix using `as.one.based()`.

  soma_dim_0_one_based <- 1 + as.numeric(tbl$GetColumnByName("soma_dim_0"))
  soma_dim_1_one_based <- 1 + as.numeric(tbl$GetColumnByName("soma_dim_1"))
  
  soma_data <- as.numeric(tbl$GetColumnByName("soma_data"))
  
  if (is.null(dims_one_based)) {
      dims_one_based <- c(max(soma_dim_0_one_based), max(soma_dim_1_one_based))
  }
  
  if(any(dims_one_based > .Machine$integer.max)) {
      error("The dimensions of the array are larger than supported by Matrix::sparseMatrix")
  }
  
  mat <- Matrix::sparseMatrix(i = soma_dim_0_one_based,
                              j = soma_dim_1_one_based,
                              x = soma_data,
                              dims = dims_one_based, repr = repr)
  matrixZeroBasedView(mat)
}

#' @description Converts a \link[Arrow]{Table} of sparse format (columns: "soma_dim_0", 
#' "soma_dim_1", "soma_data") to a \link{matrixZeroBasedView}
#' @param tbl  \link[Arrow]{Table} with columns "soma_dim_0", "soma_dim_1", and "soma_datah"
#' @param repr Optional one-character code for sparse matrix representation type
#' @param dims_one_based Numerical vectors with two elements, one for each dimension. If
#' \code{NULL}, then the following is used \code{c(max(tbl["soma_dim_0"]), max(tbl["soma_dim_1"]))}
#' @return \link{matrixZeroBasedView}
arrow_table_to_dense <- function(tbl, byrow) {

  # To instantiate the one-based Matrix::sparseMatrix, we need to add 1 to the
  # zero-based soma_dim_0 and soma_dim_1 (done by arrow_table_to_sparse). But, because these dimensions are
  # usually populated with soma_joinid, users will need to access the matrix
  # using the original, possibly-zero IDs. Therefore, we'll wrap the one-based
  # sparseMatrix with a shim providing basic access with zero-based indexes.
  # If needed, user can then explicitly ask the shim for the underlying
  # sparseMatrix using `as.one.based()`.

  nrows <- length(unique(as.numeric(tbl$GetColumnByName("soma_dim_0"))))
  ncols <- length(unique(as.numeric(tbl$GetColumnByName("soma_dim_1"))))
  soma_data <- as.numeric(tbl$GetColumnByName("soma_data"))
  
  mat <- matrix(soma_data, nrow = nrows, ncol = ncols, byrow = byrow)
  matrixZeroBasedView(mat)
}

