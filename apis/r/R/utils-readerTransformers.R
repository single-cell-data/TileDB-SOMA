soma_array_to_arrow <- function(x) {
    arrow::as_arrow_table(arrow::RecordBatch$import_from_c(x[[1]], x[[2]]))
}

arrow_table_to_sparse <- function(tbl, repr, dims = NULL) {
    
  soma_dim_0_one_based <- 1 + as.numeric(tbl$GetColumnByName("soma_dim_0"))
  soma_dim_1_one_based <- 1 + as.numeric(tbl$GetColumnByName("soma_dim_1"))
  
  soma_data <- as.numeric(tbl$GetColumnByName("soma_data"))
  
  if (is.null(dims)) {
      dims <- c(max(soma_dim_0_one_based), max(soma_dim_1_one_based))
  }
  
  if(any(dims > .Machine$integer.max)) {
      error("The dimensions of the array are larger than supported by Matrix::sparseMatrix")
  }
  
  Matrix::sparseMatrix(i = soma_dim_0_one_based,
                       j = soma_dim_1_one_based,
                       x = soma_data,
                       dims = dims, repr = repr)
  
}

