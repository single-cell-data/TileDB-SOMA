test_that("Ephemeral Experiment mechanics", {
  # Create the experiment
  uri <- withr::local_tempdir("ephemeral-experiment")
  expect_warning(EphemeralExperiment$new(uri))
  expect_no_condition(experiment <- EphemeralExperiment$new())
  expect_false(experiment$exists())
  expect_error(experiment$obs)
  expect_s3_class(experiment$create(), experiment$class())

  # Add obs
  expect_error(experiment$obs)
  obs <- create_and_populate_obs(file.path(uri, "obs"))
  expect_no_condition(experiment$obs <- obs)
  expect_equal(experiment$length(), 1)
  expect_s3_class(experiment$obs, 'SOMADataFrame')

  # Add ms
  expect_error(experiment$ms <- obs)
  expect_error(experiment$ms <- SOMAMeasurementCreate(file.path(uri, '_ms')))
  expect_no_condition(experiment$ms <- SOMACollectionCreate(file.path(uri, 'ms')))
  expect_equal(experiment$length(), 2)
  expect_s3_class(experiment$ms, 'SOMACollection')
})
