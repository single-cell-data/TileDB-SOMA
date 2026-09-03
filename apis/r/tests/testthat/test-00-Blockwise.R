test_that("Blockwise iterator for arrow tables", {
  skip_if(!extended_tests() || covr_tests())
  skip_if_not_installed("pbmc3k.tiledb") # a Suggests: pre-package 3k PBMC data
  # see https://ghrr.github.io/drat/

  tdir <- tempfile()
  tgzfile <- system.file(
    "raw-data",
    "soco-pbmc3k.tar.gz",
    package = "pbmc3k.tiledb"
  )
  untar(tarfile = tgzfile, exdir = tdir)

  uri <- file.path(tdir, "soco", "pbmc3k_processed")
  expect_true(dir.exists(uri))

  ax <- 0
  sz <- 1000L
  expqry <- SOMAExperimentOpen(uri)
  axqry <- expqry$axis_query("RNA")
  xrqry <- axqry$X("data")

  expect_error(xrqry$blockwise(axis = 2))
  expect_error(xrqry$blockwise(size = -100))

  expect_s3_class(
    bi <- xrqry$blockwise(axis = ax, size = sz, reindex_disable_on_axis = TRUE),
    "SOMASparseNDArrayBlockwiseRead"
  )

  expect_s3_class(it <- bi$tables(), "BlockwiseTableReadIter")
  expect_false(it$read_complete())

  for (i in seq.int(
    1L,
    ceiling(it$coords_axis$length() / it$coords_axis$stride)
  )) {
    at <- it$read_next()
    expect_s3_class(at, "ArrowTabular")
  }
  expect_true(it$read_complete())

  rm(bi, it, xrqry, axqry)
  axqry <- expqry$axis_query("RNA")
  xrqry <- axqry$X("data")
  bi <- xrqry$blockwise(axis = ax, size = sz, reindex_disable_on_axis = TRUE)
  it <- bi$tables()
  at <- it$concat()
  expect_s3_class(at, "Table")
  expect_s3_class(at, "ArrowTabular")
  expect_equal(dim(at), c(4848644, 3))
})

test_that("Table blockwise iterator: re-indexed", {
  skip_if(!extended_tests() || covr_tests())
  skip_if_not_installed(
    "SeuratObject",
    minimum_version = .MINIMUM_SEURAT_VERSION("c")
  )

  obj <- get_data("pbmc_small", package = "SeuratObject")
  obj <- suppressWarnings(SeuratObject::UpdateSeuratObject(obj))
  for (lyr in setdiff(SeuratObject::Layers(obj), "data")) {
    SeuratObject::LayerData(obj, lyr) <- NULL
  }
  for (reduc in SeuratObject::Reductions(obj)) {
    obj[[reduc]] <- NULL
  }
  for (grph in SeuratObject::Graphs(obj)) {
    obj[[grph]] <- NULL
  }
  for (cmd in SeuratObject::Command(obj)) {
    obj[[cmd]] <- NULL
  }

  tmp <- tempfile("blockwise-reindexed-tables")
  uri <- write_soma(obj, uri = tmp)

  exp <- SOMAExperimentOpen(uri)
  on.exit(exp$close(), add = TRUE, after = FALSE)

  ax <- 0L
  sz <- 23L
  query <- exp$axis_query("RNA")
  xrqry <- query$X("data")

  expect_s3_class(
    bi <- xrqry$blockwise(axis = ax, size = sz, reindex_disable_on_axis = NA),
    "SOMASparseNDArrayBlockwiseRead"
  )

  expect_s3_class(it <- bi$tables(), "BlockwiseTableReadIter")
  expect_false(it$read_complete())
  expect_true(it$reindexable)
  expect_error(it$concat(), class = "notConcatenatableError")
  expect_length(it$axes_to_reindex, 0L)

  for (i in seq.int(
    1L,
    ceiling(it$coords_axis$length() / it$coords_axis$stride)
  )) {
    at <- it$read_next()
    expect_true(R6::is.R6(at))
    expect_s3_class(at, "Table")
    sd0 <- at$GetColumnByName("soma_dim_0")$as_vector()
    expect_true(min(sd0) >= 0L)
    expect_true(max(sd0) <= sz)
    strider <- attr(at, "coords")$soma_dim_0
    expect_s3_class(strider, "CoordsStrider")
    expect_true(strider$start == sz * (i - 1L))
    expect_true(strider$end < sz * i)
  }

  expect_s3_class(
    bi <- suppressWarnings(xrqry$blockwise(
      axis = ax,
      size = sz,
      reindex_disable_on_axis = FALSE
    )),
    "SOMASparseNDArrayBlockwiseRead"
  )
  expect_s3_class(it <- bi$tables(), "BlockwiseTableReadIter")
  expect_false(it$read_complete())
  expect_true(it$reindexable)
  expect_error(it$concat(), class = "notConcatenatableError")
  expect_length(it$axes_to_reindex, it$array$ndim() - 1L)

  for (i in seq.int(
    1L,
    ceiling(it$coords_axis$length() / it$coords_axis$stride)
  )) {
    at <- it$read_next()
    expect_true(R6::is.R6(at))
    expect_s3_class(at, "Table")
    sd0 <- at$GetColumnByName("soma_dim_0")$as_vector()
    expect_true(min(sd0) >= 0L)
    expect_true(max(sd0) <= sz)
  }
})

test_that("Blockwise iterator for sparse matrices", {
  skip_if(!extended_tests() || covr_tests())
  skip_if_not_installed("pbmc3k.tiledb") # a Suggests: pre-package 3k PBMC data
  # see https://ghrr.github.io/drat/

  tdir <- tempfile()
  tgzfile <- system.file(
    "raw-data",
    "soco-pbmc3k.tar.gz",
    package = "pbmc3k.tiledb"
  )
  untar(tarfile = tgzfile, exdir = tdir)

  uri <- file.path(tdir, "soco", "pbmc3k_processed")
  expect_true(dir.exists(uri))

  ax <- 0
  sz <- 1000L
  expqry <- SOMAExperimentOpen(uri)
  axqry <- expqry$axis_query("RNA")
  xrqry <- axqry$X("data")

  expect_error(xrqry$blockwise(axis = 2))
  expect_error(xrqry$blockwise(size = -100))

  expect_s3_class(
    bi <- xrqry$blockwise(axis = ax, size = sz, reindex_disable_on_axis = TRUE),
    "SOMASparseNDArrayBlockwiseRead"
  )

  expect_error(bi$sparse_matrix("C"))
  expect_error(bi$sparse_matrix("R"))

  expect_s3_class(it <- bi$sparse_matrix(), "BlockwiseSparseReadIter")
  expect_false(it$reindexable)
  expect_false(it$read_complete())

  for (i in seq.int(
    1L,
    ceiling(it$coords_axis$length() / it$coords_axis$stride)
  )) {
    at <- it$read_next()
    expect_s4_class(at, "dgTMatrix")
  }
  expect_true(it$read_complete())

  rm(bi, it, xrqry, axqry)
  axqry <- expqry$axis_query("RNA")
  xrqry <- axqry$X("data")
  bi <- xrqry$blockwise(axis = ax, size = sz, reindex_disable_on_axis = TRUE)
  it <- bi$sparse_matrix()
  at <- it$concat()
  expect_s4_class(at, "dgTMatrix")
  # Re-indexing is disabled on both axes, so blocks keep global coordinates
  # and span the full array extent (for this legacy array, its domain)
  expect_equal(dim(at), as.integer(it$array$shape()))
})

test_that("Sparse matrix blockwise iterator: re-indexed", {
  skip_if(!extended_tests() || covr_tests())
  skip_if_not_installed(
    "SeuratObject",
    minimum_version = .MINIMUM_SEURAT_VERSION("c")
  )

  obj <- get_data("pbmc_small", package = "SeuratObject")
  obj <- suppressWarnings(SeuratObject::UpdateSeuratObject(obj))
  for (lyr in setdiff(SeuratObject::Layers(obj), "data")) {
    SeuratObject::LayerData(obj, lyr) <- NULL
  }
  for (reduc in SeuratObject::Reductions(obj)) {
    obj[[reduc]] <- NULL
  }
  for (grph in SeuratObject::Graphs(obj)) {
    obj[[grph]] <- NULL
  }
  for (cmd in SeuratObject::Command(obj)) {
    obj[[cmd]] <- NULL
  }

  tmp <- tempfile("blockwise-reindexed-sparse")
  uri <- write_soma(obj, uri = tmp)

  exp <- SOMAExperimentOpen(uri)
  on.exit(exp$close(), add = TRUE, after = FALSE)

  ax <- 0L
  sz <- 23L
  query <- exp$axis_query("RNA")
  xrqry <- query$X("data")

  # Reindex only on major axis
  expect_s3_class(
    bi <- xrqry$blockwise(axis = ax, size = sz, reindex_disable_on_axis = NA),
    "SOMASparseNDArrayBlockwiseRead"
  )

  expect_error(bi$sparse_matrix("C"))

  expect_s3_class(it <- bi$sparse_matrix(), "BlockwiseSparseReadIter")
  expect_false(it$read_complete())
  expect_true(it$reindexable)
  expect_error(it$concat(), class = "notConcatenatableError")
  expect_length(it$axes_to_reindex, 0L)

  for (i in seq.int(
    1L,
    ceiling(it$coords_axis$length() / it$coords_axis$stride)
  )) {
    mat <- it$read_next()
    expect_s4_class(mat, "TsparseMatrix")
    dims <- c(
      ifelse(it$read_complete(), yes = ncol(obj) %% sz, no = sz),
      nrow(obj)
    )
    expect_identical(dim(mat), dims)
    expect_true(min(mat@i) >= 0L)
    expect_true(max(mat@i) <= sz)
    strider <- attr(mat, "coords")$soma_dim_0
    expect_s3_class(strider, "CoordsStrider")
    expect_true(strider$start == sz * (i - 1L))
    expect_true(strider$end < sz * i)
  }

  # Reindex on all axes
  expect_s3_class(
    bi <- suppressWarnings(xrqry$blockwise(
      axis = ax,
      size = sz,
      reindex_disable_on_axis = FALSE
    )),
    "SOMASparseNDArrayBlockwiseRead"
  )
  expect_s3_class(it <- bi$sparse_matrix(), "BlockwiseSparseReadIter")
  expect_false(it$read_complete())
  expect_true(it$reindexable)
  expect_error(it$concat(), class = "notConcatenatableError")
  expect_length(it$axes_to_reindex, it$array$ndim() - 1L)

  for (i in seq.int(
    1L,
    ceiling(it$coords_axis$length() / it$coords_axis$stride)
  )) {
    mat <- it$read_next()
    expect_s4_class(mat, "TsparseMatrix")
    dims <- c(
      ifelse(it$read_complete(), yes = ncol(obj) %% sz, no = sz),
      nrow(obj)
    )
    expect_identical(dim(mat), dims)
    expect_true(min(mat@i) >= 0L)
    expect_true(max(mat@i) <= sz)
  }
})

test_that("Blockwise iterate through full array", {
  skip_if(!extended_tests() || covr_tests())

  uri <- tempfile("blockwise-complete")
  n_obs <- 500L
  n_var <- 210L
  X_layer <- "counts"
  exp <- create_and_populate_experiment(
    uri,
    n_obs = n_obs,
    n_var = n_var,
    X_layer_names = X_layer,
    mode = "READ"
  )

  on.exit(exp$close(), add = TRUE, after = FALSE)

  n_chunks <- 8L
  # Stride across `obs`
  obs_stride <- n_obs %/% n_chunks
  it <- exp$ms$get("RNA")$X$get(X_layer)$read()$blockwise(
    axis = 0L,
    size = obs_stride
  )$sparse_matrix()
  expect_false(it$read_complete())
  i <- 1L
  while (!it$read_complete()) {
    i <- i + 1L
    mat <- it$read_next()
    expect_s4_class(mat, "dgTMatrix")
    expect_identical(ncol(mat), n_var)
    nobs <- ifelse(
      it$read_complete(),
      yes = n_obs %% obs_stride,
      no = obs_stride
    )
    expect_identical(nrow(mat), nobs, expected.label = nobs)
    ncoords <- ifelse(
      it$read_complete(),
      yes = n_obs %% n_chunks,
      no = obs_stride
    )
    expect_identical(
      length(attr(mat, "coords")$soma_dim_0),
      ncoords,
      expected.label = ncoords
    )
  }
  expect_true(it$read_complete())

  # Stride across `var`
  var_stride <- n_var %/% n_chunks
  it <- exp$ms$get("RNA")$X$get(X_layer)$read()$blockwise(
    axis = 1L,
    size = var_stride
  )$sparse_matrix()
  expect_false(it$read_complete())
  while (!it$read_complete()) {
    mat <- it$read_next()
    expect_s4_class(mat, "dgTMatrix")
    nvar <- ifelse(
      it$read_complete(),
      yes = n_var %% var_stride,
      no = var_stride
    )
    expect_identical(ncol(mat), nvar, expected.label = nvar)
    expect_identical(nrow(mat), n_obs)
    ncoords <- ifelse(
      it$read_complete(),
      yes = n_var %% n_chunks,
      no = var_stride
    )
    expect_identical(
      length(attr(mat, "coords")$soma_dim_1),
      ncoords,
      expected.label = ncoords
    )
  }
  expect_true(it$read_complete())
})

test_that("CoordsStrider length and chunking for ranges not starting at zero", {
  strider <- CoordsStrider$new(start = 5L, end = 14L, stride = 3L)
  expect_equal(strider$length(), 10)
  expect_length(as.list(strider), 4L)
  expect_equal(
    tiledbsoma:::strider_coords(CoordsStrider$new(start = 5L, end = 14L)),
    bit64::as.integer64(5:14)
  )
  expect_true(tiledbsoma:::strider_is_full_domain(
    CoordsStrider$new(start = 0L, end = 19L),
    extent = 20L
  ))
  expect_false(tiledbsoma:::strider_is_full_domain(
    CoordsStrider$new(start = 0L, end = 18L),
    extent = 20L
  ))
  expect_false(tiledbsoma:::strider_is_full_domain(
    CoordsStrider$new(bit64::as.integer64(0:19)),
    extent = 20L
  ))
})

# Helpers for the coordinate-restricted blockwise tests below (SOMA-990).
#
# All coordinates in these helpers are 0-based SOMA coordinates. The one place
# we convert to R's 1-based indexing is inside `expected_block()`.

# A fully-populated 20x20 array where every cell holds a unique value, so a
# block that lands in the wrong place or with the wrong shape cannot pass by
# accident. Returns the URI and the dense reference matrix.
blockwise_fixture <- function() {
  n <- 20L
  mat <- Matrix::Matrix(matrix(seq_len(n * n), n, n), sparse = TRUE)
  uri <- tempfile("blockwise-restricted")
  arr <- SOMASparseNDArrayCreate(uri, arrow::int32(), shape = c(n, n))
  arr$write(mat)
  arr$close()
  list(uri = uri, mat = as.matrix(mat))
}

# Split a coordinate vector into consecutive blocks of `size`
chunk <- function(x, size) {
  unname(split(x, ceiling(seq_along(x) / size)))
}

# Read every block from a blockwise iterator into a list. `repr = NULL`
# yields Arrow tables, otherwise sparse matrices in that representation
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

# The dense matrix a block covering `rows` x `cols` of `mat` should produce.
#
# An axis that is re-indexed ("compacted") has its coordinates mapped to
# 0:(n - 1) and an extent of n, where n is the number of requested coordinates.
# An axis that is not re-indexed keeps its global coordinates and the full
# array extent. This is the rule the Python API implements.
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

# Same idea for a `tables()` block: the COO triples we expect to read back,
# sorted by coordinate
expected_table <- function(mat, rows, cols, compact_rows, compact_cols) {
  dense <- expected_block(mat, rows, cols, compact_rows, compact_cols)
  ij <- which(dense != 0, arr.ind = TRUE)
  as_coo(data.frame(
    soma_dim_0 = ij[, "row"] - 1L,
    soma_dim_1 = ij[, "col"] - 1L,
    soma_data = dense[ij]
  ))
}

# Normalize an Arrow table or data frame of COO triples for comparison
as_coo <- function(x) {
  df <- as.data.frame(x)[, c("soma_dim_0", "soma_dim_1", "soma_data")]
  df[] <- lapply(df, as.numeric)
  df <- df[order(df$soma_dim_0, df$soma_dim_1), ]
  rownames(df) <- NULL
  df
}

dense <- function(mat) unname(as.matrix(mat))

# Compare every block from `read_blocks()` against its expected value. Sparse
# matrix blocks are compared as dense matrices, table blocks as COO triples
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

# Which axes get compacted under each `reindex_disable_on_axis` setting when
# iterating over axis 0
reindex_modes <- list(
  list(reindex = NA, rows = TRUE, cols = FALSE), # default: iterated axis only
  list(reindex = FALSE, rows = TRUE, cols = TRUE), # all axes
  list(reindex = TRUE, rows = FALSE, cols = FALSE), # no axes
  list(reindex = 0L, rows = FALSE, cols = TRUE) # minor axis only
)

test_that("Blockwise tables honor a restriction on the non-iterated axis", {
  fx <- blockwise_fixture()
  arr <- SOMASparseNDArrayOpen(fx$uri)
  on.exit(arr$close(), add = TRUE, after = FALSE)
  rows <- 5:14
  cols <- 3:9
  coords <- list(soma_dim_0 = rows, soma_dim_1 = cols)
  row_chunks <- chunk(rows, 3L)

  for (mode in reindex_modes) {
    blocks <- read_blocks(arr, coords, axis = 0L, size = 3L, mode$reindex)
    expected <- lapply(row_chunks, function(br) {
      expected_table(fx$mat, br, cols, mode$rows, mode$cols)
    })
    expect_blocks(
      blocks,
      expected,
      sprintf("reindex = %s", format(mode$reindex))
    )
  }
})

test_that("Blockwise sparse matrices honor a restriction on the non-iterated axis", {
  fx <- blockwise_fixture()
  arr <- SOMASparseNDArrayOpen(fx$uri)
  on.exit(arr$close(), add = TRUE, after = FALSE)
  rows <- 5:14
  cols <- 3:9
  coords <- list(soma_dim_0 = rows, soma_dim_1 = cols)
  row_chunks <- chunk(rows, 3L)

  for (mode in reindex_modes) {
    blocks <- read_blocks(
      arr,
      coords,
      axis = 0L,
      size = 3L,
      mode$reindex,
      repr = "T"
    )
    expected <- lapply(row_chunks, function(br) {
      expected_block(fx$mat, br, cols, mode$rows, mode$cols)
    })
    expect_blocks(
      blocks,
      expected,
      sprintf("reindex = %s", format(mode$reindex))
    )
  }

  # The compressed representation from the original report
  blocks <- read_blocks(
    arr,
    coords,
    axis = 0L,
    size = 3L,
    reindex = NA,
    repr = "R"
  )
  expect_s4_class(blocks[[1L]], "dgRMatrix")
  expected <- lapply(row_chunks, function(br) {
    expected_block(fx$mat, br, cols, compact_rows = TRUE, compact_cols = FALSE)
  })
  expect_blocks(blocks, expected, "repr = R")

  # Iterate over the second axis with a restriction on the first
  col_chunks <- chunk(cols, 4L)
  blocks <- read_blocks(
    arr,
    coords,
    axis = 1L,
    size = 4L,
    reindex = NA,
    repr = "T"
  )
  expected <- lapply(col_chunks, function(bc) {
    expected_block(fx$mat, rows, bc, compact_rows = FALSE, compact_cols = TRUE)
  })
  expect_blocks(blocks, expected, "axis = 1")
})

test_that("Blockwise reads accept coords for a subset of dimensions", {
  fx <- blockwise_fixture()
  arr <- SOMASparseNDArrayOpen(fx$uri)
  on.exit(arr$close(), add = TRUE, after = FALSE)
  all_rows <- seq_len(nrow(fx$mat)) - 1L
  all_cols <- all_rows

  # Restrict only the minor axis; the iterated axis spans the full domain
  cols <- 3:9
  row_chunks <- chunk(all_rows, 8L)
  blocks <- read_blocks(
    arr,
    list(soma_dim_1 = cols),
    axis = 0L,
    size = 8L,
    reindex = NA,
    repr = "T"
  )
  expected <- lapply(row_chunks, function(br) {
    expected_block(fx$mat, br, cols, compact_rows = TRUE, compact_cols = FALSE)
  })
  expect_blocks(blocks, expected, "minor axis only")

  # Restrict only the iterated axis
  rows <- 5:14
  row_chunks <- chunk(rows, 4L)
  blocks <- read_blocks(
    arr,
    list(soma_dim_0 = rows),
    axis = 0L,
    size = 4L,
    reindex = NA,
    repr = "T"
  )
  expected <- lapply(row_chunks, function(br) {
    expected_block(
      fx$mat,
      br,
      all_cols,
      compact_rows = TRUE,
      compact_cols = FALSE
    )
  })
  expect_blocks(blocks, expected, "major axis only")
})
