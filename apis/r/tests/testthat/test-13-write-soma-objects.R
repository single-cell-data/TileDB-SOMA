test_that("write_soma.data.frame mechanics", {
  skip_if(!extended_tests())

  uri <- tempfile(pattern = "write-soma-data-frame")
  collection <- SOMACollectionCreate(uri)

  co2 <- get_data("CO2", package = "datasets")
  expect_no_condition(sdf <- write_soma(co2, uri = "co2", soma = collection))
  expect_s3_class(sdf, "SOMADataFrame")
  expect_true(sdf$exists())
  expect_identical(sdf$uri, file.path(collection$uri, "co2"))
  expect_identical(sdf$dimnames(), "soma_joinid")
  expect_identical(sdf$attrnames(), c(names(co2), "obs_id"))
  expect_error(sdf$shape(), class = "notYetImplementedError")
  schema <- sdf$schema()
  expect_s3_class(schema, "Schema")
  expect_equal(schema$num_fields - 2L, ncol(co2))
  expect_identical(
    setdiff(schema$names, c("soma_joinid", "obs_id")),
    names(co2)
  )

  collection$close()
})

test_that("write_soma.data.frame enumerations", {
  skip_if(!extended_tests())

  uri <- tempfile(pattern = "write-soma-data-frame-enumerations")
  collection <- SOMACollectionCreate(uri)
  on.exit(collection$close(), add = TRUE)

  co2 <- get_data("CO2", package = "datasets")
  expect_no_condition(sdf <- write_soma(co2, uri = "co2", soma = collection))
  expect_s3_class(sdf, "SOMADataFrame")
  expect_true(sdf$exists())
  expect_identical(sdf$uri, file.path(collection$uri, "co2"))
  sdf$close()
  expect_no_condition(sdf <- SOMADataFrameOpen(sdf$uri))

  expect_s3_class(schema <- sdf$schema(), "Schema")
  expect_equal(schema$num_fields - 2L, ncol(co2))
  expect_identical(
    setdiff(schema$names, c("soma_joinid", "obs_id")),
    names(co2)
  )

  expect_s3_class(tbl <- sdf$read()$concat(), "Table")
  expect_equal(ncol(tbl), schema$num_fields)
  expect_identical(sort(names(tbl)), sort(schema$names))

  expect_s3_class(df <- as.data.frame(tbl), "data.frame")
  for (i in names(co2)) {
    if (is.factor(co2[[i]])) {
      expect_true(is.factor(df[[i]]))
    } else {
      expect_equal(class(df[[i]]), class(co2[[i]]))
    }
  }
})

test_that("write_soma.data.frame no enumerations", {
  skip_if(!extended_tests())

  uri <- tempfile(pattern = "write-soma-data-frame-factorless")
  collection <- SOMACollectionCreate(uri)
  on.exit(collection$close(), add = TRUE)

  co2 <- get_data("CO2", package = "datasets")
  for (i in names(co2)) {
    if (is.factor(co2[[i]])) {
      co2[[i]] <- as.vector(co2[[i]])
    }
  }

  expect_no_condition(sdf <- write_soma(co2, uri = "co2", soma = collection))
  expect_s3_class(sdf, "SOMADataFrame")
  expect_true(sdf$exists())
  expect_identical(sdf$uri, file.path(collection$uri, "co2"))
  sdf$close()
  expect_no_condition(sdf <- SOMADataFrameOpen(sdf$uri))

  expect_s3_class(schema <- sdf$schema(), "Schema")
  expect_equal(schema$num_fields - 2L, ncol(co2))
  expect_identical(
    setdiff(schema$names, c("soma_joinid", "obs_id")),
    names(co2)
  )

  expect_s3_class(tbl <- sdf$read()$concat(), "Table")
  expect_equal(ncol(tbl), schema$num_fields)
  expect_identical(sort(names(tbl)), sort(schema$names))

  expect_s3_class(df <- as.data.frame(tbl), "data.frame")
  for (i in names(co2)) {
    if (is.factor(co2[[i]])) {
      expect_false(is.factor(df[[i]]))
      expect_type(df[[i]], "character")
    } else {
      expect_equal(class(df[[i]]), class(co2[[i]]))
    }
  }
})

test_that("write_soma.data.frame registration", {
  skip_if(!extended_tests())

  uri <- tempfile(pattern = "write-soma-data-frame-registration")
  collection <- SOMACollectionCreate(uri)
  on.exit(collection$close(), add = TRUE)

  co2 <- get_data("CO2", package = "datasets")

  expect_no_condition(
    sdf <- write_soma(co2, uri = "co2", soma_parent = collection, key = "co2")
  )
  expect_s3_class(sdf, "SOMADataFrame")
  expect_true(sdf$exists())
  expect_identical(sdf$uri, file.path(collection$uri, "co2"))

  sdf$close()
  collection$close()

  expect_no_condition(collection <- SOMACollectionOpen(collection$uri))
  expect_s3_class(collection, "SOMACollection")
  expect_identical(collection$length(), 1L)
  expect_identical(collection$names(), "co2")
  expect_s3_class(cdf <- collection$get("co2"), "SOMADataFrame")
  expect_s3_class(df <- as.data.frame(cdf$read()$concat()), "data.frame")
  for (col in names(co2)) {
    expect_identical(df[[col]], co2[[col]])
  }

  # Registration assertions
  expect_error(write_soma(co2, "uri", soma_parent = collection, key = TRUE))
  expect_error(write_soma(co2, "uri", soma_parent = collection, key = 1L))
  expect_error(write_soma(
    co2,
    "uri",
    soma_parent = collection,
    key = c("a", "b")
  ))
  expect_error(write_soma(co2, "uri", soma_parent = NULL, key = "co2"))
})

test_that("write_soma dense matrix mechanics", {
  skip_if(!extended_tests())

  uri <- tempfile(pattern = "write-soma-dense-matrix")
  collection <- SOMACollectionCreate(uri)

  state77 <- get(x = "state.x77", envir = getNamespace("datasets"))
  expect_no_condition(
    dmat <- write_soma(
      state77,
      uri = "state77",
      soma = collection,
      sparse = FALSE
    )
  )
  expect_s3_class(dmat, "SOMADenseNDArray")
  expect_true(dmat$exists())
  expect_identical(dmat$uri, file.path(collection$uri, "state77"))
  expect_equal(dmat$ndim(), 2L)
  expect_identical(dmat$dimnames(), paste0("soma_dim_", c(0L, 1L)))
  expect_identical(dmat$attrnames(), "soma_data")
  expect_equal(dmat$shape(), dim(state77))
  # Test transposition
  expect_no_condition(
    tmat <- write_soma(
      state77,
      uri = "state77t",
      soma = collection,
      sparse = FALSE,
      transpose = TRUE
    )
  )
  expect_s3_class(tmat, "SOMADenseNDArray")
  expect_true(tmat$exists())
  expect_identical(tmat$uri, file.path(collection$uri, "state77t"))
  expect_equal(tmat$ndim(), 2L)
  expect_identical(tmat$dimnames(), paste0("soma_dim_", c(0L, 1L)))
  expect_identical(tmat$attrnames(), "soma_data")
  expect_equal(tmat$shape(), rev(dim(state77)))
  # Error if given sparse matrix and ask for dense
  knex <- get_data("KNex", package = "Matrix")$mm
  expect_error(write_soma(
    knex,
    uri = "knex",
    soma = collection,
    sparse = FALSE
  ))
  # Work on dgeMatrices
  expect_no_condition(
    emat <- write_soma(
      as(knex, "unpackedMatrix"),
      uri = "knexd",
      soma = collection,
      sparse = FALSE
    )
  )
  expect_s3_class(emat, "SOMADenseNDArray")
  expect_true(emat$exists())
  expect_identical(emat$uri, file.path(collection$uri, "knexd"))
  expect_equal(emat$ndim(), 2L)
  expect_identical(emat$dimnames(), paste0("soma_dim_", c(0L, 1L)))
  expect_identical(emat$attrnames(), "soma_data")
  expect_equal(emat$shape(), dim(knex))

  collection$close()
})

test_that("write_soma dense matrix registration", {
  skip_if(!extended_tests())

  uri <- tempfile(pattern = "write-soma-dense-matrix-registration")
  collection <- SOMACollectionCreate(uri)
  on.exit(collection$close(), add = TRUE)

  state77 <- get(x = "state.x77", envir = getNamespace("datasets"))

  expect_no_condition(
    dmat <- write_soma(
      state77,
      uri = "state77",
      soma_parent = collection,
      sparse = FALSE,
      key = "state77"
    )
  )
  expect_s3_class(dmat, "SOMADenseNDArray")
  expect_true(dmat$exists())
  expect_identical(dmat$uri, file.path(collection$uri, "state77"))

  dmat$close()
  collection$close()

  expect_no_condition(collection <- SOMACollectionOpen(collection$uri))
  expect_s3_class(collection, "SOMACollection")
  expect_identical(collection$length(), 1L)
  expect_identical(collection$names(), "state77")
  expect_s3_class(cmat <- collection$get("state77"), "SOMADenseNDArray")
  expect_type(mat <- cmat$read_dense_matrix(), "double")
  expect_true(is.matrix(mat))
  expect_identical(mat, unname(state77))

  # Registration assertions
  expect_error(write_soma(state77, "uri", soma_parent = collection, key = TRUE))
  expect_error(write_soma(state77, "uri", soma_parent = collection, key = 1L))
  expect_error(write_soma(
    state77,
    "uri",
    soma_parent = collection,
    key = c("a", "b")
  ))
  expect_error(write_soma(state77, "uri", soma_parent = NULL, key = "state77"))
})

test_that("write_soma sparse matrix mechanics", {
  skip_if(!extended_tests())
  uri <- tempfile(pattern = "write-soma-sparse-matrix")
  collection <- SOMACollectionCreate(uri)
  knex <- get_data("KNex", package = "Matrix")$mm
  expect_no_condition(smat <- write_soma(knex, uri = "knex", soma = collection))
  expect_s3_class(smat, "SOMASparseNDArray")
  expect_true(smat$exists())
  expect_identical(smat$uri, file.path(collection$uri, "knex"))
  expect_equal(smat$ndim(), 2L)
  expect_identical(smat$dimnames(), paste0("soma_dim_", c(0L, 1L)))
  expect_identical(smat$attrnames(), "soma_data")
  expect_equal(smat$shape(), dim(knex))
  # Test transposition
  expect_no_condition(
    tmat <- write_soma(knex, uri = "knext", soma = collection, transpose = TRUE)
  )
  expect_s3_class(tmat, "SOMASparseNDArray")
  expect_true(tmat$exists())
  expect_identical(tmat$uri, file.path(collection$uri, "knext"))
  expect_equal(tmat$ndim(), 2L)
  expect_identical(tmat$dimnames(), paste0("soma_dim_", c(0L, 1L)))
  expect_identical(tmat$attrnames(), "soma_data")
  expect_equal(tmat$shape(), rev(dim(knex)))
  # Try a dense matrix
  state77 <- get(x = "state.x77", envir = getNamespace("datasets"))
  expect_no_condition(
    cmat <- write_soma(state77, "state77s", soma = collection)
  )
  expect_s3_class(cmat, "SOMASparseNDArray")
  expect_true(cmat$exists())
  expect_identical(cmat$uri, file.path(collection$uri, "state77s"))
  expect_equal(cmat$ndim(), 2L)
  expect_identical(cmat$dimnames(), paste0("soma_dim_", c(0L, 1L)))
  expect_identical(cmat$attrnames(), "soma_data")
  expect_equal(cmat$shape(), dim(state77))

  collection$close()
})

test_that("write_soma sparse matrix registration", {
  skip_if(!extended_tests())

  uri <- tempfile(pattern = "write-sparse-dense-matrix-registration")
  collection <- SOMACollectionCreate(uri)
  on.exit(collection$close(), add = TRUE)

  knex <- get_data("KNex", package = "Matrix")$mm

  expect_no_condition(
    smat <- write_soma(
      knex,
      uri = "knex",
      soma_parent = collection,
      key = "knex"
    )
  )
  expect_s3_class(smat, "SOMASparseNDArray")
  expect_true(smat$exists())
  expect_identical(smat$uri, file.path(collection$uri, "knex"))

  smat$close()
  collection$close()

  expect_no_condition(collection <- SOMACollectionOpen(collection$uri))
  expect_s3_class(collection, "SOMACollection")
  expect_identical(collection$length(), 1L)
  expect_identical(collection$names(), "knex")
  expect_s3_class(cmat <- collection$get("knex"), "SOMASparseNDArray")
  expect_s4_class(mat <- cmat$read()$sparse_matrix()$concat(), "dgTMatrix")
  expect_identical(as.matrix(mat), as.matrix(unname(knex)))

  # Registration assertions
  expect_error(write_soma(knex, "uri", soma_parent = collection, key = TRUE))
  expect_error(write_soma(knex, "uri", soma_parent = collection, key = 1L))
  expect_error(write_soma(
    knex,
    "uri",
    soma_parent = collection,
    key = c("a", "b")
  ))
  expect_error(write_soma(knex, "uri", soma_parent = NULL, key = "knex"))
})

test_that("write_soma.character mechanics", {
  skip_if(!extended_tests())

  uri <- tempfile(pattern = "write-soma-character")
  collection <- SOMACollectionCreate(uri)
  on.exit(collection$close(), add = TRUE)

  mtcars <- get_data("mtcars", package = "datasets")
  cars <- row.names(mtcars)

  expect_no_condition(sdf <- write_soma(cars, "cars", soma_parent = collection))
  expect_s3_class(sdf, "SOMADataFrame")
  expect_true(sdf$exists())
  expect_length(sdf$attrnames(), 1L)
  expect_identical(sdf$attrnames(), "values")

  hint <- "soma_uns_outgest_hint"
  expect_true(hint %in% names(sdf$get_metadata()))
  expect_identical(sdf$get_metadata(hint), uns_hint("1d")[[hint]])

  sdf$close()
  expect_no_condition(sdf <- SOMADataFrameOpen(sdf$uri))
  expect_s3_class(tbl <- sdf$read()$concat(), "Table")
  expect_identical(nrow(tbl), length(cars))
  expect_identical(tbl$values$as_vector(), cars)
})

test_that("write_soma.character scalar", {
  skip_if(!extended_tests())

  uri <- tempfile(pattern = "write-soma-character-scalar")
  collection <- SOMACollectionCreate(uri)
  on.exit(collection$close(), add = TRUE)

  mtcars <- get_data("mtcars", package = "datasets")
  cars <- paste(row.names(mtcars), collapse = ",")

  expect_no_condition(sdf <- write_soma(cars, "cars", soma_parent = collection))
  expect_s3_class(sdf, "SOMADataFrame")
  expect_true(sdf$exists())
  expect_length(sdf$attrnames(), 1L)
  expect_identical(sdf$attrnames(), "values")

  hint <- "soma_uns_outgest_hint"
  expect_true(hint %in% names(sdf$get_metadata()))
  expect_identical(sdf$get_metadata(hint), uns_hint("1d")[[hint]])

  sdf$close()
  expect_no_condition(sdf <- SOMADataFrameOpen(sdf$uri))
  expect_s3_class(tbl <- sdf$read()$concat(), "Table")
  expect_identical(nrow(tbl), 1L)
  expect_identical(tbl$soma_joinid$as_vector(), 0L)
  expect_identical(sdf.cars <- tbl$values$as_vector(), cars)
  expect_length(unlist(strsplit(sdf.cars, ",")), nrow(mtcars))
})

test_that("write_soma.IterableMatrix mechanics", {
  skip_if(!extended_tests() || covr_tests())
  skip_if_not_installed("BPCells")

  uri <- tempfile(pattern = "write-soma-bpcells")
  collection <- SOMACollectionCreate(uri)
  on.exit(collection$close(), add = TRUE, after = FALSE)

  mat <- create_sparse_matrix_with_int_dims(nrows = 1610L, ncols = 560L)
  ctx <- SOMAContext$new(c(
    soma.init_buffer_bytes = Sys.getenv(
      "TILEDB_SOMA_INIT_BUFFER_BYTES",
      unset = "33554432"
    )
  ))
  formats <- c("memory", "10x", "anndata", "dir", "hdf5")
  for (fmt in formats) {
    bpmat <- write_bpcells(mat, dirname = uri, format = fmt)
    expect_no_condition(smat <- write_soma(
      bpmat,
      uri = fmt,
      soma_parent = collection,
      context = ctx
    ))
    expect_s3_class(smat, "SOMASparseNDArray")
    expect_true(smat$exists(), info = fmt)
    expect_identical(smat$uri, file.path(collection$uri, fmt), info = fmt)
    expect_equal(smat$ndim(), 2L, info = fmt)
    expect_identical(smat$dimnames(), sprintf("soma_dim_%i", 0:1), info = fmt)
    expect_identical(smat$attrnames(), "soma_data", info = fmt)
    expect_equal(smat$shape(), dim(bpmat), info = fmt)
  }
  # Test transposition
  for (fmt in formats) {
    bpmat <- write_bpcells(mat, dirname = uri, format = fmt)
    tfmt <- sprintf("%s-transposed", fmt)
    expect_no_condition(smat <- write_soma(
      bpmat,
      uri = tfmt,
      soma_parent = collection,
      transpose = TRUE,
      context = ctx
    ))
    expect_s3_class(smat, "SOMASparseNDArray")
    expect_true(smat$exists(), info = tfmt)
    expect_identical(smat$uri, file.path(collection$uri, tfmt), info = tfmt)
    expect_equal(smat$ndim(), 2L, info = tfmt)
    expect_identical(smat$dimnames(), sprintf("soma_dim_%i", 0:1), info = tfmt)
    expect_identical(smat$attrnames(), "soma_data", info = tfmt)
    expect_equal(smat$shape(), rev(dim(bpmat)), info = tfmt)
  }
  # Test chunking
  ctx <- SOMAContext$new(c(
    soma.init_buffer_bytes = as.character(2L * (1024L ^ 2L))
  ))
  for (fmt in formats) {
    bpmat <- write_bpcells(mat, dirname = uri, format = fmt)
    cfmt <- sprintf("%s-chunked", fmt)
    expect_no_condition(smat <- write_soma(
      bpmat,
      uri = cfmt,
      soma_parent = collection,
      context = ctx
    ))
    expect_s3_class(smat, "SOMASparseNDArray")
    expect_true(smat$exists(), info = cfmt)
    expect_identical(smat$uri, file.path(collection$uri, cfmt), info = cfmt)
    expect_equal(smat$ndim(), 2L, info = cfmt)
    expect_identical(smat$dimnames(), sprintf("soma_dim_%i", 0:1), info = cfmt)
    expect_identical(smat$attrnames(), "soma_data", info = cfmt)
    expect_equal(smat$shape(), dim(bpmat), info = cfmt)
  }
})

test_that("write_soma.IterableMatrix registration", {
  skip_if(!extended_tests())
  skip_if_not_installed("BPCells")

  uri <- tempfile(pattern = "write-soma-bpcells-registration")
  collection <- SOMACollectionCreate(uri)
  on.exit(collection$close(), add = TRUE)

  mat <- create_sparse_matrix_with_int_dims(nrows = 1610L, ncols = 560L)
  ctx <- SOMAContext$new(c(
    soma.init_buffer_bytes = Sys.getenv(
      "TILEDB_SOMA_INIT_BUFFER_BYTES",
      unset = "33554432"
    )
  ))

  formats <- c("memory", "10x", "anndata", "dir", "hdf5")
  for (i in seq_along(formats)) {
    collection$reopen("WRITE")
    fmt <- formats[i]
    info <- sprintf("registration: %s", fmt)
    bpmat <- write_bpcells(mat, dirname = uri, format = fmt)
    expect_no_condition(smat <- write_soma(
      bpmat,
      uri = fmt,
      soma_parent = collection,
      key = fmt,
      context = ctx
    ))
    expect_s3_class(smat, "SOMASparseNDArray")
    expect_true(smat$exists(), info = info)
    expect_identical(smat$uri, file.path(collection$uri, fmt), info = info)

    smat$close()
    collection$reopen("READ")

    expect_s3_class(collection, "SOMACollection")
    expect_identical(collection$length(), i, info = info)
    expect_identical(collection$names(), formats[1:i], info = info)
    expect_s3_class(cmat <- collection$get(fmt), "SOMASparseNDArray")
    expect_s4_class(mat <- cmat$read()$sparse_matrix()$concat(), "dgTMatrix")
    expect_identical(as.matrix(mat), as.matrix(unname(mat)))

    # Registration assertions
    # only need to do these once
    if (i != 1L) {
      next
    }
    expect_error(write_soma(bpmat, "uri", soma_parent = collection, key = TRUE))
    expect_error(write_soma(bpmat, "uri", soma_parent = collection, key = 1L))
    expect_error(write_soma(
      bpmat,
      "uri",
      soma_parent = collection,
      key = c("a", "b")
    ))
    expect_error(write_soma(bpmat, "uri", soma_parent = NULL, key = "knex"))
  }
})

test_that("write_soma.IterableMatrix integrity", {
  skip_if(!extended_tests())
  skip_if_not_installed("BPCells")
  uri <- tempfile(pattern = "write-soma-bpcells-integrity")
  collection <- SOMACollectionCreate(uri)
  on.exit(collection$close(), add = TRUE, after = FALSE)

  mat <- create_sparse_matrix_with_int_dims(nrows = 1610L, ncols = 560L)
  ctx <- SOMAContext$new(c(
    soma.init_buffer_bytes = Sys.getenv(
      "TILEDB_SOMA_INIT_BUFFER_BYTES",
      unset = "33554432"
    )
  ))
  formats <- c("memory", "10x", "anndata", "dir", "hdf5")
  for (fmt in formats) {
    info <- sprintf("integrity: %s", fmt)
    bpmat <- write_bpcells(mat, dirname = uri, format = fmt)
    stopifnot(all.equal(
      target = suppressMessages(suppressWarnings(as.matrix(bpmat))),
      current = as.matrix(mat),
      check.attributes = FALSE
    ))
    expect_no_condition(smat <- write_soma(
      bpmat,
      uri = fmt,
      soma_parent = collection,
      key = fmt,
      context = ctx
    ))
    expect_s3_class(smat, "SOMASparseNDArray")
    expect_true(smat$exists(), info = info)
    expect_equal(smat$shape(), dim(bpmat), info = info)
    smat$reopen("READ")
    expect_s4_class(rmat <- smat$read()$sparse_matrix()$concat(), "dgTMatrix")
    expect_identical(dim(rmat), dim(mat), info = info)
    expect_identical(as.matrix(rmat), as.matrix(mat), info = info)
  }

  # Empty chunks
  stride <- .block_size(n = ncol(mat), context = ctx)
  nchunks <- ceiling(x = nrow(x = mat) / stride)
  cases <- list(
    trailing = rbind(
      mat,
      create_sparse_matrix_with_int_dims(
        nrows = (stride * nchunks) - nrow(x = mat),
        ncols = ncol(x = mat)
      ),
      Matrix::Matrix(
        data = 0L,
        nrow = stride %/% 3L,
        ncol = ncol(x = mat),
        sparse = TRUE
      )
    ),
    middle = rbind(
      mat[seq.int(from = 1L, to = stride), , drop = FALSE],
      Matrix::Matrix(0L, nrow = stride, ncol = ncol(x = mat), sparse = TRUE),
      mat[seq.int(from = stride + 1L, to = nrow(x = mat)), , drop = FALSE]
    ),
    leading = rbind(
      Matrix::Matrix(0L, nrow = stride, ncol = ncol(x = mat), sparse = TRUE),
      mat
    )
  )
  for (i in seq_along(along.with = cases)) {
    mm <- cases[[i]]
    for (fmt in formats) {
      info <- sprintf(fmt = "%s empty: %s", names(x = cases)[i], fmt)
      key <- sprintf(fmt = "%s-%s", names(x = cases)[i], fmt)
      bpmat <- write_bpcells(mm, dirname = uri, format = fmt)
      stopifnot(all.equal(
        target = suppressMessages(suppressWarnings(as.matrix(bpmat))),
        current = as.matrix(mm),
        check.attributes = FALSE
      ))
      expect_no_condition(smat <- write_soma(
        bpmat,
        uri = key,
        soma_parent = collection,
        key = key,
        context = ctx
      ))
      expect_s3_class(smat, "SOMASparseNDArray")
      expect_true(smat$exists(), info = info)
      expect_equal(smat$shape(), dim(bpmat), info = info)
      smat$reopen("READ")
      expect_s4_class(rmat <- smat$read()$sparse_matrix()$concat(), "dgTMatrix")
      expect_identical(dim(rmat), dim(mm), info = info)
      expect_identical(as.matrix(rmat), as.matrix(mm), info = info)
    }
  }
})

test_that("get_{some,tiledb}_object_type", {
  skip_if_not_installed(
    "SeuratObject",
    minimum_version = .MINIMUM_SEURAT_VERSION("c")
  )

  suppressMessages({
    library(SeuratObject)
    library(tiledbsoma)
  })

  ## write out a SOMA
  data("pbmc_small")
  uri <- tempfile()
  expect_equal(write_soma(pbmc_small, uri = uri), uri) # uri return is success

  # SOMA
  soma_context_handle = create_soma_context()
  expect_equal(
    tiledbsoma:::get_soma_object_type(uri, soma_context_handle),
    "SOMAExperiment"
  )
  expect_equal(
    tiledbsoma:::get_soma_object_type(file.path(uri, "ms/RNA"), soma_context_handle),
    "SOMAMeasurement"
  )
  coll <- c("ms", "ms/RNA/obsm", "ms/RNA/obsp/", "ms/RNA/varm")
  for (co in coll) {
    expect_equal(
      tiledbsoma:::get_soma_object_type(file.path(uri, co), soma_context_handle),
      "SOMACollection"
    )
  }
  expect_equal(
    tiledbsoma:::get_soma_object_type(
      file.path(uri, "ms/RNA/var"),
      soma_context_handle
    ),
    "SOMADataFrame"
  )
  sparr <- c("ms/RNA/obsm/X_pca", "ms/RNA/obsm/X_tsne", "ms/RNA/obsp/RNA_snn")
  for (a in sparr) {
    expect_equal(
      tiledbsoma:::get_soma_object_type(file.path(uri, a), soma_context_handle),
      "SOMASparseNDArray"
    )
  }
  expect_error(tiledbsoma:::get_some_object_type("doesnotexit", soma_context_handle))

  ## TileDB
  grps <- c("", "ms", "ms/RNA", "ms/RNA/obsm", "ms/RNA/obsp/", "ms/RNA/varm")
  for (g in grps) {
    expect_equal(
      tiledbsoma:::get_tiledb_object_type(file.path(uri, g), soma_context_handle),
      "GROUP"
    )
  }
  arrs <- c(
    "ms/RNA/obsm/X_pca",
    "ms/RNA/obsm/X_tsne",
    "ms/RNA/obsp/RNA_snn",
    "ms/RNA/var"
  )
  for (a in arrs) {
    expect_equal(
      tiledbsoma:::get_tiledb_object_type(file.path(uri, a), soma_context_handle),
      "ARRAY"
    )
  }
  expect_equal(
    tiledbsoma:::get_tiledb_object_type("doesnotexit", soma_context_handle),
    "INVALID"
  )
})

test_that("write_soma throws existingKeyWarning for duplicate keys", {
  uri <- withr::local_tempdir("write-soma-duplicate-key")

  # Test data
  df1 <- data.frame(a = 1:5, b = letters[1:5])
  df2 <- data.frame(x = 6:10, y = letters[6:10])

  # Create collection with first object
  collection <- SOMACollectionCreate(uri)
  withr::defer(collection$close())

  sdf1 <- write_soma(df1, uri = "foo", soma_parent = collection, key = "foo")
  sdf1$close()
  expect_true("foo" %in% collection$names())

  # Attempt to write another object with same key
  # With verbose = TRUE, warning should be thrown
  expect_warning(
    withr::with_options(
      list(verbose = TRUE),
      sdf2 <- write_soma(df2, uri = "bar", soma_parent = collection, key = "foo")
    ),
    class = "existingKeyWarning"
  )

  # Verify original data is still there (not replaced)
  expect_true("foo" %in% collection$names())
  expect_equal(collection$length(), 1L)

  collection$close()
  collection <- SOMACollectionOpen(uri)
  expect_identical(
    collection$get("foo")$read()$concat()$a$as_vector(),
    df1$a
  )
})
