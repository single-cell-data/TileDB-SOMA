
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

  co2 <- get_data('CO2', package = 'datasets')
})

test_that("Resume-mode sparse arrays", {
  skip_if(!extended_tests())

  mat <- get_data('KNex', package = 'Matrix')$mm
})

test_that("Resume-mode dense arrays", {
  skip_if(!extended_tests())
  skip_if_not_installed('datasets')

  mat <- get(x = 'state.x77', envir = getNamespace('datasets'))
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
