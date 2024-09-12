
test_that("Blockwise iterator for arrow tables", {
  skip_if(!extended_tests() || covr_tests())
  skip_if_not_installed("pbmc3k.tiledb")      # a Suggests: pre-package 3k PBMC data
  # see https://ghrr.github.io/drat/

  tdir <- tempfile()
  tgzfile <- system.file("raw-data", "soco-pbmc3k.tar.gz", package="pbmc3k.tiledb")
  untar(tarfile = tgzfile, exdir = tdir)

  uri <- file.path(tdir, "soco", "pbmc3k_processed")
  expect_true(dir.exists(uri))

  ax <- 0
  sz <- 1000L
  expqry <- SOMAExperimentOpen(uri)
  axqry <- expqry$axis_query("RNA")
  xrqry <- axqry$X("data")

  expect_error(xrqry$blockwise(axis=2))
  expect_error(xrqry$blockwise(size=-100))

  expect_s3_class(
    bi <- xrqry$blockwise(axis=ax, size=sz, reindex_disable_on_axis = TRUE),
    "SOMASparseNDArrayBlockwiseRead"
  )

  expect_s3_class(it <- bi$tables(), "BlockwiseTableReadIter")
  expect_false(it$read_complete())

  for (i in seq.int(1L, ceiling(it$coords_axis$length() / it$coords_axis$stride))) {
    at <- it$read_next()
    expect_s3_class(at, "ArrowTabular")
  }
  expect_true(it$read_complete())

  rm(bi, it, xrqry, axqry)
  axqry <- expqry$axis_query("RNA")
  xrqry <- axqry$X("data")
  bi <- xrqry$blockwise(axis=ax, size=sz, reindex_disable_on_axis = TRUE)
  it <- bi$tables()
  at <- it$concat()
  expect_s3_class(at, "Table")
  expect_s3_class(at, "ArrowTabular")
  expect_equal(dim(at), c(4848644, 3))
})

test_that("Table blockwise iterator: re-indexed", {
  skip_if(!extended_tests() || covr_tests())
  skip_if_not_installed("SeuratObject", minimum_version = .MINIMUM_SEURAT_VERSION('c'))

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

  for (i in seq.int(1L, ceiling(it$coords_axis$length() / it$coords_axis$stride))) {
    at <- it$read_next()
    expect_true(R6::is.R6(at))
    expect_s3_class(at, "Table")
    sd0 <- at$GetColumnByName("soma_dim_0")$as_vector()
    expect_true(min(sd0) >= 0L)
    expect_true(max(sd0) <= sz)
    strider <- attr(at, 'coords')$soma_dim_0
    expect_s3_class(strider, 'CoordsStrider')
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

  for (i in seq.int(1L, ceiling(it$coords_axis$length() / it$coords_axis$stride))) {
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
  skip_if_not_installed("pbmc3k.tiledb")      # a Suggests: pre-package 3k PBMC data
  # see https://ghrr.github.io/drat/

  tdir <- tempfile()
  tgzfile <- system.file("raw-data", "soco-pbmc3k.tar.gz", package="pbmc3k.tiledb")
  untar(tarfile = tgzfile, exdir = tdir)

  uri <- file.path(tdir, "soco", "pbmc3k_processed")
  expect_true(dir.exists(uri))

  ax <- 0
  sz <- 1000L
  expqry <- SOMAExperimentOpen(uri)
  axqry <- expqry$axis_query("RNA")
  xrqry <- axqry$X("data")

  expect_error(xrqry$blockwise(axis=2))
  expect_error(xrqry$blockwise(size=-100))

  expect_s3_class(
    bi <- xrqry$blockwise(axis=ax, size=sz, reindex_disable_on_axis = TRUE),
    "SOMASparseNDArrayBlockwiseRead"
  )

  expect_error(bi$sparse_matrix("C"))
  expect_error(bi$sparse_matrix("R"))

  expect_s3_class(it <- bi$sparse_matrix(), "BlockwiseSparseReadIter")
  expect_false(it$read_complete())

  for (i in seq.int(1L, ceiling(it$coords_axis$length() / it$coords_axis$stride))) {
    at <- it$read_next()
    expect_s4_class(at, "dgTMatrix")
  }
  expect_true(it$read_complete())

  rm(bi, it, xrqry, axqry)
  axqry <- expqry$axis_query("RNA")
  xrqry <- axqry$X("data")
  bi <- xrqry$blockwise(axis=ax, size=sz, reindex_disable_on_axis = TRUE)
  it <- bi$sparse_matrix()
  at <- it$concat()
  expect_s4_class(at, "dgTMatrix")
  expect_equal(dim(at), c(2638, 1838))
})

test_that("Sparse matrix blockwise iterator: re-indexed", {
  skip_if(!extended_tests() || covr_tests())
  skip_if_not_installed("SeuratObject", minimum_version = .MINIMUM_SEURAT_VERSION('c'))

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

  for (i in seq.int(1L, ceiling(it$coords_axis$length() / it$coords_axis$stride))) {
    mat <- it$read_next()
    expect_s4_class(mat, "TsparseMatrix")
    expect_identical(dim(mat), rev(dim(obj)))
    expect_true(min(mat@i) >= 0L)
    expect_true(max(mat@i) <= sz)
    strider <- attr(mat, 'coords')$soma_dim_0
    expect_s3_class(strider, 'CoordsStrider')
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
  expect_s3_class(it <- bi$sparse_matrix(), "BlockwiseSparseReadIter")
  expect_false(it$read_complete())
  expect_true(it$reindexable)
  expect_error(it$concat(), class = "notConcatenatableError")
  expect_length(it$axes_to_reindex, it$array$ndim() - 1L)

  for (i in seq.int(1L, ceiling(it$coords_axis$length() / it$coords_axis$stride))) {
    mat <- it$read_next()
    expect_s4_class(mat, "TsparseMatrix")
    expect_identical(dim(mat), rev(dim(obj)))
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
  it <- exp$ms$get("RNA")$X$get(X_layer)$read()$blockwise(axis = 0L, size = obs_stride)$sparse_matrix()
  expect_false(it$read_complete())
  while (!it$read_complete()) {
    mat <- it$read_next()
    expect_s4_class(mat, "dgTMatrix")
    expect_identical(ncol(mat), n_var)
    expect_identical(nrow(mat), n_obs)
    expect_identical(
      length(attr(mat, "coords")$soma_dim_0),
      ifelse(it$read_complete(), yes = n_obs %% n_chunks, no = obs_stride)
    )
  }
  expect_true(it$read_complete())

  # Stride across `var`
  var_stride <- n_var %/% n_chunks
  it <- exp$ms$get("RNA")$X$get(X_layer)$read()$blockwise(axis = 1L, size = var_stride)$sparse_matrix()
  expect_false(it$read_complete())
  while (!it$read_complete()) {
    mat <- it$read_next()
    expect_s4_class(mat, "dgTMatrix")
    expect_identical(ncol(mat), n_var)
    expect_identical(nrow(mat), n_obs)
    expect_identical(
      length(attr(mat, "coords")$soma_dim_1),
      ifelse(it$read_complete(), yes = n_var %% n_chunks, no = var_stride)
    )
  }
  expect_true(it$read_complete())
})
