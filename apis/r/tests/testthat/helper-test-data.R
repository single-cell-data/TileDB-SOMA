create_sparse_matrix_with_int_dims <- function(nrows = 10, ncols = 5, seed = 1, repr = "T") {
  set.seed(seed)
  Matrix::rsparsematrix(
    nrow = nrows,
    ncol = ncols,
    density = 0.6,
    rand.x = function(n) as.integer(runif(n, min = 1, max = 100)),
    repr = repr
  )
}

create_sparse_matrix_with_string_dims <- function(nrows = 10, ncols = 5, seed = 1, repr = "T") {
  smat <- create_sparse_matrix_with_int_dims(nrows, ncols, seed, repr)
  dimnames(smat) <- list(
    paste0("i", seq_len(nrows)),
    paste0("j", seq_len(ncols))
  )
  smat
}

create_dense_matrix_with_int_dims <- function(nrows = 10, ncols = 5, seed = 1) {
  set.seed(seed)
  matrix(
    data = as.integer(runif(nrows * ncols, min = 1, max = 100)),
    nrow = nrows,
    ncol = ncols
  )
}

create_arrow_schema <- function() {
  arrow::schema(
    arrow::field("foo", arrow::int32(), nullable = FALSE),
    arrow::field("soma_joinid", arrow::int64(), nullable = FALSE),
    arrow::field("bar", arrow::float64(), nullable = FALSE),
    arrow::field("baz", arrow::large_utf8(), nullable = FALSE)
  )
}

create_and_populate_soma_dataframe <- function(uri) {

  arrow_schema <- create_arrow_schema()

  tbl <- arrow::arrow_table(
    foo = 1L:10L,
    soma_joinid = 1L:10L,
    bar = 1.1:10.1,
    baz = letters[1:10],
    schema = arrow_schema
  )

  sdf <- SOMADataFrame$new(uri)
  sdf$create(arrow_schema, index_column_names = "foo")
  sdf$write(tbl)
  sdf
}
