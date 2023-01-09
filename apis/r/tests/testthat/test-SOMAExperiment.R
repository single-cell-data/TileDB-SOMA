test_that("Basic mechanics", {
  uri <- withr::local_tempdir("soma-experiment")


  experiment <- SOMAExperiment$new(uri)

  expect_false(experiment$exists())
  expect_error(experiment$obs, "Group does not exist.")

  experiment$create()

  # Add obs
  expect_error(experiment$obs, "No member named 'obs' found")
  obs <- create_and_populate_obs(file.path(uri, "obs"))
  experiment$obs <- obs
  expect_true(inherits(experiment$obs, "SOMADataFrame"))
})
