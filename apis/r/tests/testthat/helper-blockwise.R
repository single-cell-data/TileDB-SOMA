#' Split a coordinate vector into consecutive blocks
#'
#' @param x A vector of coordinates
#' @param size The number of coordinates per block
#'
#' @return A list of vectors, each up to \code{size} elements of \code{x},
#' split in consecutive order
#'
chunk <- function(x, size) {
  unname(split(x, ceiling(seq_along(x) / size)))
}

#' Create a densely-filled reference \code{SOMASparseNDArray} fixture for the
#' coordinate-restricted blockwise tests (SOMA-990)
#'
#' @return A list with:
#' \itemize{
#'   \item \code{uri}: the URI of a 20x20 \code{SOMASparseNDArray}, densely
#'   filled so every cell holds a unique value, so a block that lands in the
#'   wrong place or with the wrong shape cannot pass by accident
#'   \item \code{mat}: the dense reference matrix written to the array
#' }
#'
#' @note All coordinates used with this fixture are 0-based SOMA
#' coordinates; the one place we convert to R's 1-based indexing is inside
#' \code{expected_block()}
#'
blockwise_fixture <- function() {
  n <- 20L
  mat <- Matrix::Matrix(matrix(seq_len(n * n), n, n), sparse = TRUE)
  uri <- tempfile("blockwise-restricted")
  arr <- SOMASparseNDArrayCreate(uri, arrow::int32(), shape = c(n, n))
  arr$write(mat)
  arr$close()
  list(uri = uri, mat = as.matrix(mat))
}

#' Read every block from a blockwise iterator into a list
#'
#' @param arr An open \code{SOMASparseNDArray}
#' @param coords Coordinates to pass to \code{arr$read()}
#' @param axis The axis to iterate over blockwise
#' @param size The block size
#' @param reindex The \code{reindex_disable_on_axis} argument to
#' \code{blockwise()}
#' @param repr If \code{NULL}, yields Arrow tables via \code{tables()};
#' otherwise the sparse matrix representation to yield via
#' \code{sparse_matrix()}
#'
#' @return A list of blocks read from the iterator, in order
#'
read_blocks <- function(arr, coords, axis, size, reindex, repr = NULL) {
  bi <- arr$read(coords)$blockwise(
    axis = axis,
    size = size,
    reindex_disable_on_axis = reindex
  )
  it <- if (is.null(repr)) bi$tables() else bi$sparse_matrix(repr)
  blocks <- list()
  while (!it$read_complete()) {
    blocks[[length(blocks) + 1L]] <- it$read_next()
  }
  blocks
}

#' The dense matrix a block covering \code{rows} x \code{cols} of \code{mat}
#' should produce
#'
#' @param mat The dense reference matrix (see \code{blockwise_fixture()})
#' @param rows,cols 0-based row/column coordinates covered by this block
#' @param compact_rows,compact_cols Whether the row/column axis is
#' re-indexed ("compacted"); a re-indexed axis has its coordinates mapped to
#' \code{0:(n - 1)} and an extent of \code{n}, where \code{n} is the number
#' of requested coordinates, while a non-re-indexed axis keeps its global
#' coordinates and the full array extent. This is the rule the Python API
#' implements
#'
#' @return The dense matrix expected for this block
#'
expected_block <- function(mat, rows, cols, compact_rows, compact_cols) {
  n <- nrow(mat)
  i <- if (compact_rows) seq_along(rows) - 1L else rows
  j <- if (compact_cols) seq_along(cols) - 1L else cols
  out <- matrix(
    0,
    if (compact_rows) length(rows) else n,
    if (compact_cols) length(cols) else n
  )
  out[i + 1L, j + 1L] <- mat[rows + 1L, cols + 1L]
  out
}

#' The COO triples a \code{tables()} block covering \code{rows} x
#' \code{cols} of \code{mat} should produce
#'
#' @inheritParams expected_block
#'
#' @return A data frame of \code{soma_dim_0}, \code{soma_dim_1},
#' \code{soma_data} triples, sorted by coordinate (see \code{as_coo()})
#'
expected_table <- function(mat, rows, cols, compact_rows, compact_cols) {
  dense <- expected_block(mat, rows, cols, compact_rows, compact_cols)
  ij <- which(dense != 0, arr.ind = TRUE)
  as_coo(data.frame(
    soma_dim_0 = ij[, "row"] - 1L,
    soma_dim_1 = ij[, "col"] - 1L,
    soma_data = dense[ij]
  ))
}

#' Normalize an Arrow table or data frame of COO triples for comparison
#'
#' @param x An Arrow table or data frame with \code{soma_dim_0},
#' \code{soma_dim_1}, and \code{soma_data} columns
#'
#' @return A data frame with those three columns, coerced to numeric and
#' sorted by \code{soma_dim_0} then \code{soma_dim_1}
#'
as_coo <- function(x) {
  df <- as.data.frame(x)[, c("soma_dim_0", "soma_dim_1", "soma_data")]
  df[] <- lapply(df, as.numeric)
  df <- df[order(df$soma_dim_0, df$soma_dim_1), ]
  rownames(df) <- NULL
  df
}

#' Coerce a sparse matrix block to a plain dense matrix for comparison
#'
#' @param mat A \code{Matrix}-derived sparse matrix
#'
#' @return \code{mat} as an unnamed dense \code{matrix}
#'
dense <- function(mat) unname(as.matrix(mat))

#' Compare every block from \code{read_blocks()} against its expected value
#'
#' @param blocks A list of blocks, as returned by \code{read_blocks()}
#' @param expected A list of the same length as \code{blocks}, giving the
#' expected value of each block (see \code{expected_block()} /
#' \code{expected_table()})
#' @param label A label included in test failure messages to identify which
#' scenario failed
#'
#' @return \code{NULL}, invisibly; called for its \code{testthat::expect_*()}
#' side effects. Sparse matrix blocks are compared as dense matrices (see
#' \code{dense()}), table blocks as COO triples (see \code{as_coo()})
#'
expect_blocks <- function(blocks, expected, label) {
  expect_length(blocks, length(expected))
  for (i in seq_along(blocks)) {
    actual <- if (inherits(blocks[[i]], "Matrix")) {
      dense(blocks[[i]])
    } else {
      as_coo(blocks[[i]])
    }
    expect_equal(
      actual,
      expected[[i]],
      info = sprintf("%s, block %d", label, i)
    )
  }
}

#' Which axes get compacted under each \code{reindex_disable_on_axis}
#' setting, when iterating over axis 0
#'
#' @format A list of scenarios, each a list with:
#' \itemize{
#'   \item \code{reindex}: the \code{reindex_disable_on_axis} value to pass
#'   to \code{blockwise()}
#'   \item \code{rows}, \code{cols}: whether the row (iterated) / column
#'   (minor) axis is compacted under that value
#' }
#'
reindex_modes <- list(
  list(reindex = NA, rows = TRUE, cols = FALSE), # default: iterated axis only
  list(reindex = FALSE, rows = TRUE, cols = TRUE), # all axes
  list(reindex = TRUE, rows = FALSE, cols = FALSE), # no axes
  list(reindex = 0L, rows = FALSE, cols = TRUE) # minor axis only
)
