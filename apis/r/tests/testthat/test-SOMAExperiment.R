test_that("Basic mechanics", {
  uri <- withr::local_tempdir("soma-experiment")


  experiment <- SOMAExperiment$new(uri)

  expect_false(experiment$exists())
  expect_error(experiment$obs, "Group does not exist.")

  experiment$create()
  # TODO: Determine behavior for retrieving empty obs/ms
  # expect_null(experiment$obs)
  # expect_null(experiment$ms)

  # Add obs
  expect_error(experiment$obs, "No member named 'obs' found")

  obs <- create_and_populate_obs(file.path(uri, "obs"))
  experiment$obs <- obs
  expect_equal(experiment$length(), 1)
  expect_true(inherits(experiment$obs, "SOMADataFrame"))

  # Add ms
  expect_error(experiment$ms <- obs, "ms must be a 'SOMAMeasurement'")

  measurement <- SOMAMeasurement$new(file.path(uri, "ms"))
  measurement$create()
  experiment$ms <- measurement
  expect_equal(experiment$length(), 2)
  expect_true(inherits(experiment$ms, "SOMAMeasurement"))
})
