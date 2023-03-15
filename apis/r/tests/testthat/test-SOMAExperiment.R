test_that("Basic mechanics", {
  uri <- withr::local_tempdir("soma-experiment")

  experiment <- SOMAExperiment$new(uri, internal_use_only = "allowed_use")

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
  expect_error(experiment$ms <- obs, "ms must be a 'SOMACollection'")
  expect_error(
    experiment$ms <- SOMAMeasurement$new(file.path(uri, "_ms"),
                                         internal_use_only = "allowed_use")$create(),
    "ms must be a 'SOMACollection'"
  )

  experiment$ms <- SOMACollection$new(file.path(uri, "ms"),
                                      internal_use_only = "allowed_use")$create()
  expect_equal(experiment$length(), 2)
  expect_true(inherits(experiment$ms, "SOMACollection"))
})

test_that("Configured SOMAExperiment", {
  cfg <- PlatformConfig$new()
  cfg$set(
    'tiledb',
    'create',
    value = ScalarMap$new()$setv(
      capacity = 8888L,
      cell_order = 'row-major',
      tile_order = 'col-major'
    )
  )
  uri <- withr::local_tempdir("soma-experiment-config")
  n_obs <- 20L
  n_var <- 10L
  experiment <- create_and_populate_experiment(
    uri = uri,
    n_obs = n_obs,
    n_var = n_var,
    X_layer_names = c("counts", "logcounts"),
    config = cfg
  )
  expect_equal(experiment$tiledbsoma_ctx$get('capacity'), '8888')
  expect_equal(experiment$tiledbsoma_ctx$get('cell_order'), 'row-major')
  expect_equal(experiment$tiledbsoma_ctx$get('tile_order'), 'col-major')
})
