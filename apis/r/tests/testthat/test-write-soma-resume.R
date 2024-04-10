
factories <- list(
  substitute(SOMADataFrameCreate),
  substitute(SOMASparseNDArrayCreate),
  substitute(SOMADenseNDArrayCreate),
  substitute(SOMACollectionCreate),
  substitute(SOMAMeasurementCreate),
  substitute(SOMAExperimentCreate)
)

schema <- arrow::infer_schema(data.frame(
  soma_joinid = bit64::integer64(),
  int = integer()
))

test_that("Factory re-creation", {
  skip_if(!extended_tests())
  for (i in seq_along(factories)) {
    fname <- as.character(factories[[i]])
    fxn <- eval(factories[[i]])
    uri <- withr::local_tempdir(fname)
    expect_no_condition(obj <- switch(
      EXPR = fname,
      SOMADataFrameCreate = fxn(uri, schema = schema),
      SOMASparseNDArrayCreate = ,
      SOMADenseNDArrayCreate = fxn(uri, type = arrow::int32(), shape = c(20L, 10L)),
      fxn(uri)
    ))
    expect_error(fxn(uri))
    obj$close()
  }
})

test_that("Resume-mode factories", {
  skip_if(!extended_tests())
  for (i in seq_along(factories)) {
    fname <- as.character(factories[[i]])
    if (fname == 'SOMADenseNDArrayCreate') {
      next
    }
    fxn <- eval(factories[[i]])
    label <- paste0(fname, "-resume")
    uri <- withr::local_tempdir(label)
    # Do an initial create
    expect_no_condition(obj <- switch(
      EXPR = fname,
      SOMADataFrameCreate = fxn(uri, schema = schema),
      SOMASparseNDArrayCreate = ,
      SOMADenseNDArrayCreate = fxn(uri, type = arrow::int32(), shape = c(20L, 10L)),
      fxn(uri)
    ))
    expect_true(obj$is_open(), label = fname)
    expect_identical(obj$mode(), "WRITE", label = fname)
    expect_true(obj$exists(), label = fname)
    obj$close()

    # Test that re-creating in "resume" mode simply re-opens the object
    expect_no_condition(obj <- switch(
      EXPR = fname,
      SOMADataFrameCreate = fxn(uri, schema = schema, ingest_mode = "resume"),
      SOMASparseNDArrayCreate = ,
      SOMADenseNDArrayCreate = fxn(
        uri,
        type = arrow::int32(),
        shape = c(20L, 10L),
        ingest_mode = "resume"
      ),
      fxn(uri, ingest_mode = "resume")
    ))
    expect_true(obj$is_open(), label = label)
    expect_identical(obj$mode(), "WRITE", label = label)
    obj$close()
  }
})

test_that("Resume-mode data frames", {
  skip_if(!extended_tests())
  skip_if_not_installed('datasets')

  collection <- SOMACollectionCreate(withr::local_tempdir("dataframe-resume"))
  on.exit(collection$close(), add = TRUE, after = FALSE)

  co2 <- get_data('CO2', package = 'datasets')

  # Test resume-mode when writing data.frames
  uri <- "co2-complete"
  expect_s3_class(
    sdf <- write_soma(co2, uri = uri, soma_parent = collection),
    "SOMADataFrame"
  )
  on.exit(sdf$close(), add = TRUE, after = FALSE)

  sdf$reopen("READ")
  df <- as.data.frame(sdf$read()$concat())
  for (i in names(co2)) {
    expect_identical(
      df[[i]],
      co2[[i]],
      label = sprintf("df[['%s']]", i),
      expected.label = sprintf("co2[['%s']]", i)
    )
  }

  # Expect error when writing to existing array
  expect_error(write_soma(co2, uri = uri, soma_parent = collection))

  # Expect seamless pass when resuming writing to exisitng array
  expect_s3_class(
    sdfr <- write_soma(
      co2,
      uri = uri,
      soma_parent = collection,
      ingest_mode = "resume"
    ),
    "SOMADataFrame"
  )
  on.exit(sdfr$close(), add = TRUE, after = FALSE)
  expect_identical(sdf$uri, sdfr$uri)

  sdfr$reopen("READ")
  dfr <- as.data.frame(sdfr$read()$concat())
  for (i in names(co2)) {
    for (i in names(co2)) {
      expect_identical(
        dfr[[i]],
        co2[[i]],
        label = sprintf("dfr[['%s']]", i),
        expected.label = sprintf("co2[['%s']]", i)
      )
    }
  }

  # Test resume-mode with partial writes
  idx <- seq.int(1L, floor(nrow(co2) / 3))
  co2p <- co2[idx, , drop = FALSE]
  uri <- "co2-parital"
  expect_s3_class(
    sdfp <- write_soma(co2p, uri = uri, soma_parent = collection),
    "SOMADataFrame"
  )
  on.exit(sdfp$close(), add = TRUE, after = FALSE)

  sdfp$reopen("READ")
  dfp <- as.data.frame(sdfp$read()$concat())
  for (i in names(co2)) {
    for (i in names(co2)) {
      expect_identical(
        dfp[[i]],
        co2[[i]][idx],
        label = sprintf("dfp[['%s']]", i),
        expected.label = sprintf("co2[['%s']]", i)
      )
      expect_false(
        identical(dfp[[i]], co2[[i]]),
        info = "partial write is not identical to original data",
        label = sprintf("dfp[['%s']] != co2[['%s']]", i, i)
      )
    }
  }

  expect_s3_class(
    sdfc <- write_soma(
      co2,
      uri = uri,
      soma_parent = collection,
      ingest_mode = "resume"
    ),
    "SOMADataFrame"
  )
  on.exit(sdfc$close(), add = TRUE, after = FALSE)
  expect_identical(sdfc$uri, sdfp$uri)

  sdfc$reopen("READ")
  dfc <- as.data.frame(sdfc$read()$concat())
  for (i in names(co2)) {
    for (i in names(co2)) {
      expect_identical(
        dfc[[i]],
        co2[[i]],
        label = sprintf("dfc[['%s']]", i),
        expected.label = sprintf("co2[['%s']]", i)
      )
    }
  }
})

test_that("Resume-mode sparse arrays", {
  skip_if(!extended_tests())

  collection <- SOMACollectionCreate(withr::local_tempdir("sparse-array-resume"))
  on.exit(collection$close(), add = TRUE, after = FALSE)

  knex <- as(get_data('KNex', package = 'Matrix')$mm, "TsparseMatrix")

  # Test resume-mode when writing sparse arrays
  uri <- "knex-complete"
  expect_s3_class(
    ssa <- write_soma(knex, uri = uri, soma_parent = collection),
    "SOMASparseNDArray"
  )
  on.exit(ssa$close(), add = TRUE, after = FALSE)

  ssa$reopen("READ")
  mat <- ssa$read()$sparse_matrix()$concat()
  expect_equal(
    as.matrix(mat),
    as.matrix(knex),
    label = "mat",
    expected.label = "knex"
  )

  # Expect error when writing to existing array
  expect_error(write_soma(knex, uri = uri, soma_parent = collection))

  # Expect seamless pass when resuming writing to complete existing array
  expect_s3_class(
    ssar <- write_soma(
      knex,
      uri = uri,
      soma_parent = collection,
      ingest_mode = "resume"
    ),
    "SOMASparseNDArray"
  )
  on.exit(ssar$close(), add = TRUE, after = FALSE)
  expect_identical(ssar$uri, ssa$uri)

  ssar$reopen("READ")
  matr <- ssar$read()$sparse_matrix()$concat()
  expect_identical(
    as.matrix(matr),
    as.matrix(knex),
    label = "matr",
    expected.label = "knex"
  )

  # Test resume-mode with partial writes
  idx <- seq.int(1L, floor(nrow(knex) / 3))
  knexp <- knex[idx, , drop = FALSE]
  uri <- "knex-parital"
  expect_s3_class(
    ssap <- write_soma(knexp, uri = uri, soma_parent = collection, shape = dim(knex)),
    "SOMASparseNDArray"
  )
  on.exit(ssap$close(), add = TRUE, after = FALSE)

  ssap$reopen("READ")
  matp <- ssap$read()$sparse_matrix()$concat()
  expect_identical(dim(matp), dim(knex))
  expect_identical(max(matp@i) + 1L, max(idx), label = "max(matp@i)")
  expect_equal(
    as.matrix(matp[idx, , drop = FALSE]),
    as.matrix(knexp),
    label = "matp",
    expected.label = "knexp"
  )

  expect_s3_class(
    ssac <- write_soma(
      knex,
      uri = uri,
      soma_parent = collection,
      ingest_mode = "resume"
    ),
    "SOMASparseNDArray"
  )
  on.exit(ssac$close(), add = TRUE, after = FALSE)
  expect_identical(ssac$uri, ssap$uri)

  ssac$reopen("READ")
  matc <- ssac$read()$sparse_matrix()$concat()
  expect_identical(dim(matc), dim(knex))
  expect_identical(max(matc@i) + 1L, nrow(knex))
  expect_identical(
    as.matrix(matc),
    as.matrix(knex),
    label = "matc",
    expected.label = "knex"
  )
  bbox <- tryCatch(
    as.integer(ssac$used_shape(simplify = TRUE, index1 = TRUE)),
    error = function(...) NULL
  )
  if (!is.null(bbox)) {
    expect_identical(bbox, dim(knex))
  }
})

test_that("Resume-mode dense arrays", {
  skip_if(!extended_tests())
  skip_if_not_installed('datasets')

  collection <- SOMACollectionCreate(withr::local_tempdir("dense-array-resume"))
  on.exit(collection$close(), add = TRUE, after = FALSE)

  mat <- get(x = 'state.x77', envir = getNamespace('datasets'))

  # Resume mode should always fail for dense arrays
  expect_s3_class(
    sda <- write_soma(
      mat,
      uri = "state-x77",
      soma_parent = collection,
      sparse = FALSE
    ),
    "SOMADenseNDArray"
  )
  on.exit(sda$close(), add = TRUE, after = FALSE)

  expect_error(write_soma(
    mat,
    uri = "state-x77",
    soma_parent = collection,
    sparse = FALSE
  ))
  expect_error(write_soma(
    mat,
    uri = "state-x77",
    soma_parent = collection,
    sparse = FALSE,
    ingest_mode = "resume"
  ))
})

test_that("Resume-mode Seurat", {
  skip_if(!extended_tests())
  skip_if_not_installed("SeuratObject", .MINIMUM_SEURAT_VERSION("c"))

  pbmc_small <- get_data("pbmc_small", package = "SeuratObject")

  # Test resume-mode when writing Seurat object
  expect_type(
    uri <- write_soma(
      pbmc_small,
      uri = withr::local_tempdir(SeuratObject::Project(pbmc_small))
    ),
    "character"
  )
  exp <- SOMAExperimentOpen(uri)
  on.exit(exp$close(), add = TRUE, after = FALSE)

  layers <- exp$ms$get("RNA")$X$names()
  names(layers) <- gsub("_", ".", layers)

  obj <- suppressWarnings(SOMAExperimentAxisQuery$new(exp, "RNA")$to_seurat(
    layers,
    obs_index = "obs_id",
    var_index = "var_id"
  ))

  expect_identical(
    sort(names(obj)),
    sort(names(pbmc_small)),
    label = "names(obj)",
    expected.label = "names(pbmc_small)"
  )

  for (i in names(obj)) {
    expect_s4_class(obj[[i]], class(pbmc_small[[i]]))
    expect_identical(dim(obj[[i]]), dim(pbmc_small[[i]]))
  }

  for (i in names(pbmc_small[[]])) {
    expect_identical(
      obj[[i]],
      pbmc_small[[i]],
      label = sprintf("obj[['%s']]", i),
      expected.label = sprintf("pbmc_small[['%s']]", i)
    )
  }

  # Expect error when writing to existing array
  expect_error(write_soma(pbmc_small, uri = uri))

  # Expect seamless pass when resuming writing to exisitng experiment
  expect_type(
    urir <- write_soma(pbmc_small, uri = uri, ingest_mode = "resume"),
    "character"
  )
  expect_identical(urir, uri)

  expr <- SOMAExperimentOpen(urir)
  on.exit(expr$close(), add = TRUE, after = FALSE)

  objr <- suppressWarnings(SOMAExperimentAxisQuery$new(expr, "RNA")$to_seurat(
    layers,
    obs_index = "obs_id",
    var_index = "var_id"
  ))

  expect_identical(
    sort(names(objr)),
    sort(names(pbmc_small)),
    label = "names(objr)",
    expected.label = "names(pbmc_small)"
  )

  for (i in names(objr)) {
    expect_s4_class(objr[[i]], class(pbmc_small[[i]]))
    expect_identical(dim(objr[[i]]), dim(pbmc_small[[i]]))
  }

  for (i in names(pbmc_small[[]])) {
    expect_identical(
      objr[[i]],
      pbmc_small[[i]],
      label = sprintf("objr[['%s']]", i),
      expected.label = sprintf("pbmc_small[['%s']]", i)
    )
  }

  expr$close()
  gc()

  # Test resume-mode with partial writes
  idx <- seq.int(1L, floor(ncol(pbmc_small) / 3))
  pbmc_partial <- subset(pbmc_small, cells = idx)
  for (i in names(pbmc_partial)) {
    if (inherits(i, "Assay")) {
      next
    }
    SeuratObject::DefaultAssay(pbmc_partial[[i]]) <- SeuratObject::DefaultAssay(pbmc_partial[[i]]) %||%
      SeuratObject::DefaultAssay(pbmc_partial)
  }

  expect_type(
    urip <- write_soma(
      pbmc_partial,
      uri = withr::local_tempdir("pbmc-partial"),
      shape = dim(pbmc_small)
    ),
    "character"
  )

  expp <- SOMAExperimentOpen(urip)
  on.exit(expp$close(), add = TRUE, after = FALSE)

  layers <- expp$ms$get("RNA")$X$names()
  names(layers) <- gsub("_", ".", layers)

  objp <- suppressWarnings(SOMAExperimentAxisQuery$new(expp, "RNA")$to_seurat(
    layers,
    obs_index = "obs_id",
    var_index = "var_id"
  ))

  expect_identical(dim(objp), dim(pbmc_partial))

  for (i in names(objp)) {
    expect_s4_class(objp[[i]], class(pbmc_partial[[i]]))
    expect_identical(dim(objp[[i]]), dim(pbmc_partial[[i]]))
  }

  for (i in names(pbmc_partial[[]])) {
    expect_identical(
      objp[[i]],
      pbmc_partial[[i]],
      label = sprintf("objp[['%s']]", i),
      expected.label = sprintf("pbmc_partial[['%s']]", i)
    )
  }

  expp$close()
  gc()

  expect_type(
    uric <- write_soma(pbmc_small, uri = urip, ingest_mode = "resume"),
    "character"
  )
  expect_identical(uric, urip)

  expc <- SOMAExperimentOpen(uric)
  on.exit(expc$close(), add = TRUE, after = FALSE)

  objc <- suppressWarnings(SOMAExperimentAxisQuery$new(expc, "RNA")$to_seurat(
    layers,
    obs_index = "obs_id",
    var_index = "var_id"
  ))

  expect_identical(dim(objc), dim(pbmc_small))

  for (i in names(objc)) {
    expect_s4_class(objc[[i]], class(pbmc_small[[i]]))
    expect_identical(dim(objc[[i]]), dim(pbmc_small[[i]]))
  }

  for (i in names(pbmc_small[[]])) {
    expect_identical(
      objc[[i]],
      pbmc_small[[i]],
      label = sprintf("objc[['%s']]", i),
      expected.label = sprintf("pbmc_small[['%s']]", i)
    )
  }
})

test_that("Resume-mode SingleCellExperiment", {
  skip_if(!extended_tests())
  skip_if_not_installed("pbmc3k.sce")
  suppressMessages(skip_if_not_installed("SingleCellExperiment", .MINIMUM_SCE_VERSION("c")))

  sce <- get_data('pbmc3k.final', package = "pbmc3k.sce")
  SingleCellExperiment::mainExpName(sce) <- "RNA"

  # Test resume-mode when writing Seurat object
  expect_type(
    uri <- write_soma(sce, uri = withr::local_tempdir("single-cell-experiment")),
    "character"
  )
  exp <- SOMAExperimentOpen(uri)
  on.exit(exp$close(), add = TRUE, after = FALSE)

  obj <- SOMAExperimentAxisQuery$new(exp, "RNA")$to_single_cell_experiment(
    obs_index = "obs_id",
    var_index = "var_id"
  )

  expect_identical(dim(obj), dim(sce))
  expect_identical(
    sort(SummarizedExperiment::assayNames(obj)),
    sort(SummarizedExperiment::assayNames(sce)),
    label = 'assayNames(obj)',
    expected.label = 'assayNames(sce)'
  )
  for (i in SummarizedExperiment::assayNames(sce)) {
    expect_identical(
      dim(SummarizedExperiment::assay(obj, i)),
      dim(SummarizedExperiment::assay(sce, i)),
      label = sprintf("dim(assay(obj, '%s'))", i),
      expected.label = sprintf("dim(assay(sce, '%s'))", i)
    )
  }

  expect_identical(
    sort(SingleCellExperiment::reducedDimNames(obj)),
    sort(SingleCellExperiment::reducedDimNames(sce)),
    label = 'reducedDimNames(obj)',
    expected.label = 'reducedDimNames(sce)'
  )
  for (i in SingleCellExperiment::reducedDimNames(sce)) {
    expect_identical(
      dim(SingleCellExperiment::reducedDim(obj, i)),
      dim(SingleCellExperiment::reducedDim(sce, i)),
      label = sprintf("dim(reducedDim(obj, '%s'))", i),
      expected.label = sprintf("dim(reducedDim(sce, '%s'))", i)
    )
  }

  expect_identical(
    sort(SingleCellExperiment::colPairNames(obj)),
    sort(SingleCellExperiment::colPairNames(sce)),
    label = 'colPairNames(obj)',
    expected.label = 'colPairNames(sce)'
  )
  for (i in SingleCellExperiment::colPairNames(sce)) {
    expect_identical(
      dim(SingleCellExperiment::colPair(obj, i)),
      dim(SingleCellExperiment::colPair(sce, i)),
      label = sprintf("dim(colPair(obj, '%s'))", i),
      expected.label = sprintf("dim(colPair(sce, '%s'))", i)
    )
  }

  expect_identical(
    sort(SingleCellExperiment::rowPairNames(obj)),
    sort(SingleCellExperiment::rowPairNames(sce)),
    label = 'rowPairNames(obj)',
    expected.label = 'rowPairNames(sce)'
  )
  for (i in SingleCellExperiment::rowPairNames(sce)) {
    expect_identical(
      dim(SingleCellExperiment::rowPair(obj, i)),
      dim(SingleCellExperiment::rowPair(sce, i)),
      label = sprintf("dim(rowPair(obj, '%s'))", i),
      expected.label = sprintf("dim(rowPair(sce, '%s'))", i)
    )
  }

  for (i in names(SummarizedExperiment::colData(sce))) {
    expect_equivalent(
      SummarizedExperiment::colData(obj)[[i]],
      SummarizedExperiment::colData(sce)[[i]],
      label = sprintf("colData(obj)[['%s']]", i),
      expected.label = sprintf("colData(sce)[['%s']]", i)
    )
  }

  for (i in names(SummarizedExperiment::rowData(sce))) {
    expect_equivalent(
      SummarizedExperiment::rowData(obj)[[i]],
      SummarizedExperiment::rowData(sce)[[i]],
      label = sprintf("rowData(obj)[['%s']]", i),
      expected.label = sprintf("rowData(sce)[['%s']]", i)
    )
  }

  exp$close()
  gc()

  # Expect error when writing to existing array
  expect_error(write_soma(sce, uri = uri))

  # Expect seamless pass when resuming writing to exisitng experiment
  expect_type(
    urir <- write_soma(sce, uri = uri, ingest_mode = "resume"),
    "character"
  )
  expect_identical(urir, uri)

  expr <- SOMAExperimentOpen(urir)
  on.exit(expr$close(), add = TRUE, after = FALSE)

  objr <- SOMAExperimentAxisQuery$new(expr, "RNA")$to_single_cell_experiment(
    obs_index = "obs_id",
    var_index = "var_id"
  )

  expect_identical(
    sort(SummarizedExperiment::assayNames(objr)),
    sort(SummarizedExperiment::assayNames(sce)),
    label = 'assayNames(objr)',
    expected.label = 'assayNames(sce)'
  )
  for (i in SummarizedExperiment::assayNames(sce)) {
    expect_identical(
      dim(SummarizedExperiment::assay(objr, i)),
      dim(SummarizedExperiment::assay(sce, i)),
      label = sprintf("dim(assay(objr, '%s'))", i),
      expected.label = sprintf("dim(assay(sce, '%s'))", i)
    )
  }

  expect_identical(
    sort(SingleCellExperiment::reducedDimNames(objr)),
    sort(SingleCellExperiment::reducedDimNames(sce)),
    label = 'reducedDimNames(objr)',
    expected.label = 'reducedDimNames(sce)'
  )
  for (i in SingleCellExperiment::reducedDimNames(sce)) {
    expect_identical(
      dim(SingleCellExperiment::reducedDim(objr, i)),
      dim(SingleCellExperiment::reducedDim(sce, i)),
      label = sprintf("dim(reducedDim(objr, '%s'))", i),
      expected.label = sprintf("dim(reducedDim(sce, '%s'))", i)
    )
  }

  expect_identical(
    sort(SingleCellExperiment::colPairNames(objr)),
    sort(SingleCellExperiment::colPairNames(sce)),
    label = 'colPairNames(objr)',
    expected.label = 'colPairNames(sce)'
  )
  for (i in SingleCellExperiment::colPairNames(sce)) {
    expect_identical(
      dim(SingleCellExperiment::colPair(objr, i)),
      dim(SingleCellExperiment::colPair(sce, i)),
      label = sprintf("dim(colPair(objr, '%s'))", i),
      expected.label = sprintf("dim(colPair(sce, '%s'))", i)
    )
  }

  expect_identical(
    sort(SingleCellExperiment::rowPairNames(objr)),
    sort(SingleCellExperiment::rowPairNames(sce)),
    label = 'rowPairNames(objr)',
    expected.label = 'rowPairNames(sce)'
  )
  for (i in SingleCellExperiment::rowPairNames(sce)) {
    expect_identical(
      dim(SingleCellExperiment::rowPair(objr, i)),
      dim(SingleCellExperiment::rowPair(sce, i)),
      label = sprintf("dim(rowPair(objr, '%s'))", i),
      expected.label = sprintf("dim(rowPair(sce, '%s'))", i)
    )
  }

  for (i in names(SummarizedExperiment::colData(sce))) {
    expect_equivalent(
      SummarizedExperiment::colData(objr)[[i]],
      SummarizedExperiment::colData(sce)[[i]],
      label = sprintf("colData(objr)[['%s']]", i),
      expected.label = sprintf("colData(sce)[['%s']]", i)
    )
  }

  for (i in names(SummarizedExperiment::rowData(sce))) {
    expect_equivalent(
      SummarizedExperiment::rowData(objr)[[i]],
      SummarizedExperiment::rowData(sce)[[i]],
      label = sprintf("rowData(objr)[['%s']]", i),
      expected.label = sprintf("rowData(sce)[['%s']]", i)
    )
  }

  expr$close()
  gc()

  # Test resume-mode with partial writes
  idx <- seq.int(1L, floor(ncol(sce) / 3))
  sce_partial <- sce[, idx]

  expect_type(
    urip <- write_soma(
      sce_partial,
      uri = withr::local_tempdir("single-cell-experiment-partial"),
      shape = dim(sce)
    ),
    "character"
  )

  expp <- SOMAExperimentOpen(urip)
  on.exit(expp$close(), add = TRUE, after = FALSE)

  objp <- SOMAExperimentAxisQuery$new(expp, "RNA")$to_single_cell_experiment(
    obs_index = "obs_id",
    var_index = "var_id"
  )

  expect_identical(dim(objp), dim(sce_partial))

  expect_identical(
    sort(SummarizedExperiment::assayNames(objp)),
    sort(SummarizedExperiment::assayNames(sce_partial)),
    label = 'assayNames(objp)',
    expected.label = 'assayNames(sce_partial)'
  )
  for (i in SummarizedExperiment::assayNames(sce_partial)) {
    expect_identical(
      dim(SummarizedExperiment::assay(objp, i)),
      dim(SummarizedExperiment::assay(sce_partial, i)),
      label = sprintf("dim(assay(objp, '%s'))", i),
      expected.label = sprintf("dim(assay(sce_partial, '%s'))", i)
    )
  }

  expect_identical(
    sort(SingleCellExperiment::reducedDimNames(objp)),
    sort(SingleCellExperiment::reducedDimNames(sce_partial)),
    label = 'reducedDimNames(objp)',
    expected.label = 'reducedDimNames(sce_partial)'
  )
  for (i in SingleCellExperiment::reducedDimNames(sce_partial)) {
    expect_identical(
      dim(SingleCellExperiment::reducedDim(objp, i)),
      dim(SingleCellExperiment::reducedDim(sce_partial, i)),
      label = sprintf("dim(reducedDim(objp, '%s'))", i),
      expected.label = sprintf("dim(reducedDim(sce_partial, '%s'))", i)
    )
  }

  expect_identical(
    sort(SingleCellExperiment::colPairNames(objp)),
    sort(SingleCellExperiment::colPairNames(sce_partial)),
    label = 'colPairNames(objp)',
    expected.label = 'colPairNames(sce_partial)'
  )
  for (i in SingleCellExperiment::colPairNames(sce_partial)) {
    expect_identical(
      dim(SingleCellExperiment::colPair(objp, i)),
      dim(SingleCellExperiment::colPair(sce_partial, i)),
      label = sprintf("dim(colPair(objp, '%s'))", i),
      expected.label = sprintf("dim(colPair(sce_partial, '%s'))", i)
    )
  }

  expect_identical(
    sort(SingleCellExperiment::rowPairNames(objp)),
    sort(SingleCellExperiment::rowPairNames(sce_partial)),
    label = 'rowPairNames(objp)',
    expected.label = 'rowPairNames(sce_partial)'
  )
  for (i in SingleCellExperiment::rowPairNames(sce_partial)) {
    expect_identical(
      dim(SingleCellExperiment::rowPair(objp, i)),
      dim(SingleCellExperiment::rowPair(sce_partial, i)),
      label = sprintf("dim(rowPair(objp, '%s'))", i),
      expected.label = sprintf("dim(rowPair(sce_partial, '%s'))", i)
    )
  }

  for (i in names(SummarizedExperiment::colData(sce_partial))) {
    expect_equivalent(
      SummarizedExperiment::colData(objp)[[i]],
      SummarizedExperiment::colData(sce_partial)[[i]],
      label = sprintf("colData(objp)[['%s']]", i),
      expected.label = sprintf("colData(sce_partial)[['%s']]", i)
    )
  }

  for (i in names(SummarizedExperiment::rowData(sce_partial))) {
    expect_equivalent(
      SummarizedExperiment::rowData(objp)[[i]],
      SummarizedExperiment::rowData(sce_partial)[[i]],
      label = sprintf("rowData(objp)[['%s']]", i),
      expected.label = sprintf("rowData(sce_partial)[['%s']]", i)
    )
  }

  expp$close()
  gc()

  expect_type(
    uric <- write_soma(sce, uri = urip, ingest_mode = "resume"),
    "character"
  )
  expect_identical(uric, urip)

  expc <- SOMAExperimentOpen(uric)
  on.exit(expc$close(), add = TRUE, after = FALSE)

  objc <- SOMAExperimentAxisQuery$new(expc, "RNA")$to_single_cell_experiment(
    obs_index = "obs_id",
    var_index = "var_id"
  )

  expect_identical(dim(objc), dim(sce))

  expect_identical(
    sort(SummarizedExperiment::assayNames(objc)),
    sort(SummarizedExperiment::assayNames(sce)),
    label = 'assayNames(objc)',
    expected.label = 'assayNames(sce)'
  )
  for (i in SummarizedExperiment::assayNames(sce)) {
    expect_identical(
      dim(SummarizedExperiment::assay(objc, i)),
      dim(SummarizedExperiment::assay(sce, i)),
      label = sprintf("dim(assay(objc, '%s'))", i),
      expected.label = sprintf("dim(assay(sce, '%s'))", i)
    )
  }

  expect_identical(
    sort(SingleCellExperiment::reducedDimNames(objc)),
    sort(SingleCellExperiment::reducedDimNames(sce)),
    label = 'reducedDimNames(objc)',
    expected.label = 'reducedDimNames(sce)'
  )
  for (i in SingleCellExperiment::reducedDimNames(sce)) {
    expect_identical(
      dim(SingleCellExperiment::reducedDim(objc, i)),
      dim(SingleCellExperiment::reducedDim(sce, i)),
      label = sprintf("dim(reducedDim(objc, '%s'))", i),
      expected.label = sprintf("dim(reducedDim(sce, '%s'))", i)
    )
  }

  expect_identical(
    sort(SingleCellExperiment::colPairNames(objc)),
    sort(SingleCellExperiment::colPairNames(sce)),
    label = 'colPairNames(objc)',
    expected.label = 'colPairNames(sce)'
  )
  for (i in SingleCellExperiment::colPairNames(sce)) {
    expect_identical(
      dim(SingleCellExperiment::colPair(objc, i)),
      dim(SingleCellExperiment::colPair(sce, i)),
      label = sprintf("dim(colPair(objc, '%s'))", i),
      expected.label = sprintf("dim(colPair(sce, '%s'))", i)
    )
  }

  expect_identical(
    sort(SingleCellExperiment::rowPairNames(objc)),
    sort(SingleCellExperiment::rowPairNames(sce)),
    label = 'rowPairNames(objc)',
    expected.label = 'rowPairNames(sce)'
  )
  for (i in SingleCellExperiment::rowPairNames(sce)) {
    expect_identical(
      dim(SingleCellExperiment::rowPair(objc, i)),
      dim(SingleCellExperiment::rowPair(sce, i)),
      label = sprintf("dim(rowPair(objc, '%s'))", i),
      expected.label = sprintf("dim(rowPair(sce, '%s'))", i)
    )
  }

  for (i in names(SummarizedExperiment::colData(sce))) {
    expect_equivalent(
      SummarizedExperiment::colData(objc)[[i]],
      SummarizedExperiment::colData(sce)[[i]],
      label = sprintf("colData(objc)[['%s']]", i),
      expected.label = sprintf("colData(sce)[['%s']]", i)
    )
  }

  for (i in names(SummarizedExperiment::rowData(sce))) {
    expect_equivalent(
      SummarizedExperiment::rowData(objc)[[i]],
      SummarizedExperiment::rowData(sce)[[i]],
      label = sprintf("rowData(objc)[['%s']]", i),
      expected.label = sprintf("rowData(sce)[['%s']]", i)
    )
  }
})
