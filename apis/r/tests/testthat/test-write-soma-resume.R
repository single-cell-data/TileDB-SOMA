
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
  on.exit(sdfc$close())
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

  mat <- get_data('KNex', package = 'Matrix')$mm
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
  skip_if(TRUE)
  skip_if(!extended_tests())
  skip_if_not_installed('SeuratObject', .MINIMUM_SEURAT_VERSION('c'))

  pbmc_small <- get_data('pbmc_small', package = 'SeuratObject')

  uri <- withr::local_tempdir(SeuratObject::Project(pbmc_small))
})

test_that("Resume-mode SingleCellExperiment", {
  skip_if(TRUE)
  skip_if(!extended_tests())
  skip_if_not_installed('pbmc3k.sce')
  suppressMessages(skip_if_not_installed('SingleCellExperiment', .MINIMUM_SCE_VERSION('c')))

  sce <- get_data('pbmc3k.final', package = 'pbmc3k.sce')
  SingleCellExperiment::mainExpName(sce) <- 'RNA'

  uri <- withr::local_tempdir('single-cell-experiment')
})
