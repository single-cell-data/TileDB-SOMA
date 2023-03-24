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

create_and_populate_soma_dataframe <- function(
  uri,
  nrows = 10L,
  seed = 1,
  index_column_names = "foo"
) {
  set.seed(seed)

  arrow_schema <- create_arrow_schema()

  tbl <- arrow::arrow_table(
    foo = seq.int(nrows) + 1000L,
    soma_joinid = bit64::seq.integer64(from = 0L, to = nrows - 1L),
    bar = seq(nrows) + 0.1,
    baz = as.character(seq.int(nrows) + 1000L),
    schema = arrow_schema
  )

  sdf <- SOMADataFrame$new(uri, internal_use_only = "allowed_use")
  sdf$create(arrow_schema, index_column_names = index_column_names)
  sdf$write(tbl)
  sdf
}

#' @param ... Arguments passed to create_sparse_matrix_with_int_dims
create_and_populate_sparse_nd_array <- function(uri, ...) {

  smat <- create_sparse_matrix_with_int_dims(...)

  ndarray <- SOMASparseNDArray$new(uri, internal_use_only = "allowed_use")
  ndarray$create(arrow::int32(), shape = dim(smat))
  ndarray$write(smat)
  ndarray
}

get_data <- function(x, package = NULL) {
  stopifnot(
    is_scalar_character(x),
    is.null(package) || is_scalar_character(package)
  )
  e <- new.env()
  on.exit(rm(e), add = TRUE)
  data(list = x, package = package, envir = e)
  return(e[[x]])
}
