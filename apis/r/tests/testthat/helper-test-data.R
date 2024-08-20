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

create_arrow_schema <- function(foo_first = TRUE) {
  bl <- FALSE
  if (foo_first) {
    arrow::schema(
      arrow::field("foo", arrow::int32(), nullable = bl),
      arrow::field("soma_joinid", arrow::int64(), nullable = bl),
      arrow::field("bar", arrow::float64(), nullable = bl),
      arrow::field("baz", arrow::large_utf8(), nullable = bl)
    )
  } else {
    arrow::schema(
      arrow::field("soma_joinid", arrow::int64(), nullable = bl),
      arrow::field("foo", arrow::int32(), nullable = bl),
      arrow::field("bar", arrow::float64(), nullable = bl),
      arrow::field("baz", arrow::large_utf8(), nullable = bl)
    )
  }
}

create_arrow_table <- function(nrows = 10L, factors = FALSE) {
  if (isTRUE(factors)) {
    return(arrow::arrow_table(
      foo = seq.int(nrows) + 1000L,
      soma_joinid = bit64::seq.integer64(from = 0L, to = nrows - 1L),
      bar = seq(nrows) + 0.1,
      baz = as.character(seq.int(nrows) + 1000L),
      grp = factor(c(
        rep_len("lvl1", length.out = floor(nrows / 2)),
        rep_len("lvl2", length.out = ceiling(nrows / 2))
      ))
    ))
  }
  arrow::arrow_table(
      foo = seq.int(nrows) + 1000L,
      soma_joinid = bit64::seq.integer64(from = 0L, to = nrows - 1L),
      bar = seq(nrows) + 0.1,
      baz = as.character(seq.int(nrows) + 1000L)
      # schema = create_arrow_schema()
    )
}
