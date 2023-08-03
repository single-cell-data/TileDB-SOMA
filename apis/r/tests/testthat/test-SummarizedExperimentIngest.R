test_that("Write SummarizedExperiment mechanics", {
  skip_if_not_installed('SummarizedExperiment', '1.28.0')
  skip_if_not_installed('airway')
  se <- get_data('airway', package = 'airway')
  uri <- withr::local_tempdir('airway')

  expect_no_condition(uri <- write_soma(se, uri, 'airway'))

  expect_type(uri, 'character')
  expect_true(grepl('^airway', basename(uri)))

  expect_no_condition(experiment <- SOMAExperimentOpen(uri))
  on.exit(experiment$close())

  expect_no_error(experiment$ms)

  expect_equal(experiment$mode(), "READ")
  expect_s3_class(experiment, 'SOMAExperiment')
  expect_true(grepl('^airway', basename(experiment$uri)))

  expect_identical(experiment$ms$names(), 'airway')
  expect_s3_class(ms <- experiment$ms$get('airway'), 'SOMAMeasurement')

  expect_identical(ms$X$names(), SummarizedExperiment::assayNames(se))

  expect_error(ms$obsm)
  expect_error(ms$varm)
  expect_error(ms$obsp)
  expect_error(ms$varp)

  expect_identical(
    setdiff(experiment$obs$attrnames(), 'obs_id'),
    names(SummarizedExperiment::colData(se))
  )

  # Test ms_name assertions
  expect_error(write_soma(se, uri))
  expect_error(write_soma(se, uri, ''))
  expect_error(write_soma(se, uri, NA_character_))
  expect_error(write_soma(se, uri, c('a', 'b')))
  expect_error(write_soma(se, uri, 1))
  expect_error(write_soma(se, uri, TRUE))
})
