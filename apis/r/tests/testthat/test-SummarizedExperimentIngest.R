test_that("Write SummarizedExperiment mechanics", {
  skip_if_not_installed('SummarizedExperiment', '1.28.0')
  skip_if_not_installed('pbmc3k.sce')

  se <- get_data('pbmc3k.final', package = 'pbmc3k.sce')
  var_df <- SummarizedExperiment::rowData(se)
  features <- rownames(se)

  se <- as(se, 'SummarizedExperiment')
  SummarizedExperiment::rowData(se) <- var_df
  rownames(se) <- features

  uri <- withr::local_tempdir('summarized-experiment')

  expect_no_condition(uri <- suppressMessages(write_soma(se, uri, 'RNA')))

  expect_type(uri, 'character')
  expect_true(grepl('^summarized-experiment', basename(uri)))

  expect_no_condition(experiment <- SOMAExperimentOpen(uri))
  expect_s3_class(experiment, 'SOMAExperiment')
  on.exit(experiment$close())

  expect_no_error(experiment$ms)

  expect_equal(experiment$mode(), "READ")
  expect_s3_class(experiment, 'SOMAExperiment')
  expect_true(grepl('^summarized-experiment', basename(experiment$uri)))

  expect_s3_class(experiment$obs, 'SOMADataFrame')

  expect_identical(experiment$ms$names(), 'RNA')
  expect_s3_class(ms <- experiment$ms$get('RNA'), 'SOMAMeasurement')

  expect_s3_class(ms$var, 'SOMADataFrame')

  expect_identical(
    sort(ms$X$names()),
    sort(SummarizedExperiment::assayNames(se))
  )

  expect_error(ms$obsm)
  expect_error(ms$varm)
  expect_error(ms$obsp)
  expect_error(ms$varp)

  expect_identical(
    setdiff(experiment$obs$attrnames(), 'obs_id'),
    names(SummarizedExperiment::colData(se))
  )

  expect_identical(
    setdiff(ms$var$attrnames(), 'var_id'),
    names(SummarizedExperiment::rowData(se))
  )

  # Test ms_name assertions
  expect_error(write_soma(se, uri))
  expect_error(write_soma(se, uri, ''))
  expect_error(write_soma(se, uri, NA_character_))
  expect_error(write_soma(se, uri, c('a', 'b')))
  expect_error(write_soma(se, uri, 1))
  expect_error(write_soma(se, uri, TRUE))
})
