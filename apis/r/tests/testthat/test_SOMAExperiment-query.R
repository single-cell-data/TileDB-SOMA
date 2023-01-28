test_that("AxisQuery", {
  query <- AxisQuery$new()
  expect_null(query$value_filter)
  expect_null(query$coords)

  query <- AxisQuery$new(value_filter = "foo")
  expect_equal(query$value_filter, "foo")
  expect_null(query$coords)

  query <- AxisQuery$new(coords = list(foo = 1L:2L))
  expect_null(query$value_filter)
  expect_equal(query$coords, list(foo = 1L:2L))

  query <- AxisQuery$new(value_filter = "foo", coords = list(foo = 1L:2L))
  expect_equal(query$value_filter, "foo")
  expect_equal(query$coords, list(foo = 1L:2L))

  # Expected failures
  expect_error(
    AxisQuery$new(value_filter = 1),
    "'value_filter' must be a scalar character"
  )
  expect_error(
    AxisQuery$new(value_filter = c("foo", "bar")),
    "'value_filter' must be a scalar character"
  )
  expect_error(
    AxisQuery$new(coords = list(1L:2L)),
    "'coords' must be a named list"
  )
})

test_that("ExperimentAxisQuery returns all coordinates by default", {
  uri <- withr::local_tempdir("soma-experiment-query-all")
  n_obs <- 20L
  n_var <- 10L

  experiment <- create_and_populate_experiment(
    uri = uri,
    n_obs = n_obs,
    n_var = n_var,
    X_layer_names = "RNA"
  )

  query <- ExperimentAxisQuery$new(
    experiment = experiment,
    measurement_name = "RNA"
  )

  expect_true(is_arrow_table(query$obs()))
  expect_true(is_arrow_table(query$var()))

  expect_true(inherits(query$obs_joinids(), "Array"))
  expect_true(inherits(query$var_joinids(), "Array"))

  expect_equal(query$n_obs, n_obs)
  expect_equal(query$n_vars, n_var)

  expect_error(query$X(), "Must specify a layer name")
  expect_error(query$X("foo"), "The following layer does not exist: foo")
})
