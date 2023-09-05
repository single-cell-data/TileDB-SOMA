test_that("Write SingleCellExperiment mechanics", {
  skip_if(!extended_tests() || covr_tests())
  skip_if_not_installed('pbmc3k.sce')
  suppressMessages(skip_if_not_installed('SingleCellExperiment', .MINIMUM_SCE_VERSION('c')))

  sce <- get_data('pbmc3k.final', package = 'pbmc3k.sce')
  SingleCellExperiment::mainExpName(sce) <- 'RNA'

  uri <- withr::local_tempdir('single-cell-experiment')

  expect_no_condition(uri <- suppressMessages(write_soma(sce, uri)))

  expect_type(uri, 'character')
  expect_true(grepl('^single-cell-experiment', basename(uri)))

  expect_no_condition(experiment <- SOMAExperimentOpen(uri))
  expect_s3_class(experiment, 'SOMAExperiment')
  on.exit(experiment$close())

  expect_no_error(experiment$ms)

  expect_identical(experiment$ms$names(), SingleCellExperiment::mainExpName(sce))
  expect_s3_class(
    ms <- experiment$ms$get(SingleCellExperiment::mainExpName(sce)),
    'SOMAMeasurement'
  )

  expect_identical(
    sort(ms$X$names()),
    sort(SummarizedExperiment::assayNames(sce))
  )

  expect_no_error(ms$obsm)
  expect_identical(
    sort(ms$obsm$names()),
    sort(SingleCellExperiment::reducedDimNames(sce))
  )
  for (i in ms$obsm$names()) {
    expect_s3_class(arr <- ms$obsm$get(i), 'SOMASparseNDArray')
    expect_equivalent(arr$shape(), dim(SingleCellExperiment::reducedDim(sce, i)))
  }

  # SCE doesn't support `varm`
  expect_error(ms$varm)

  expect_no_error(ms$obsp)
  expect_identical(
    sort(ms$obsp$names()),
    sort(SingleCellExperiment::colPairNames(sce))
  )
  for (i in ms$obsp$names()) {
    expect_s3_class(arr <- ms$obsp$get(i), 'SOMASparseNDArray')
    expect_equivalent(arr$shape(), dim(SingleCellExperiment::colPair(sce, i)))
  }

  expect_no_error(ms$varp)
  expect_identical(
    sort(ms$varp$names()),
    sort(SingleCellExperiment::rowPairNames(sce))
  )
  for (i in ms$varp$names()) {
    expect_s3_class(arr <- ms$varp$get(i), 'SOMASparseNDArray')
    expect_equivalent(arr$shape(), dim(SingleCellExperiment::rowPair(sce, i)))
  }
})

test_that("SingleCellExperiment mainExpName mechanics", {
  skip_if(!extended_tests() || covr_tests())
  skip_if_not_installed('SingleCellExperiment', .MINIMUM_SCE_VERSION('c'))
  skip_if_not_installed('pbmc3k.sce')

  sce <- get_data('pbmc3k.final', package = 'pbmc3k.sce')
  SingleCellExperiment::mainExpName(sce) <- NULL
  expect_null(SingleCellExperiment::mainExpName(sce))

  uri <- withr::local_tempdir('mainexp')
  ms_name <- 'RNA'

  expect_error(uri <- suppressMessages(write_soma(sce, uri)))

  expect_no_condition(uri <- suppressMessages(write_soma(sce, uri, ms_name)))
  expect_no_condition(experiment <- SOMAExperimentOpen(uri))
  expect_identical(experiment$ms$names(), ms_name)

  uri <- withr::local_tempdir('mainexp-2')
  ms_name2 <- 'ASSAY'
  SingleCellExperiment::mainExpName(sce) <- 'MY_NAME'

  expect_type(SingleCellExperiment::mainExpName(sce), 'character')

  expect_no_condition(uri <- suppressMessages(write_soma(sce, uri, ms_name2)))
  expect_no_condition(experiment <- SOMAExperimentOpen(uri))
  expect_identical(experiment$ms$names(), ms_name2)
})
