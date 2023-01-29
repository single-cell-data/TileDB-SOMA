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
    X_layer_names = c("counts", "logcounts")
  )

  query <- ExperimentAxisQuery$new(
    experiment = experiment,
    measurement_name = "RNA"
  )

  # obs/var tables
  expect_true(query$obs()$Equals(experiment$obs$read()))
  expect_true(query$var()$Equals(experiment$ms$get("RNA")$var$read()))

  # obs/var joinids
  expect_equal(
    query$obs_joinids(),
    arrow::concat_arrays(experiment$obs$read()$soma_joinid)
  )
  expect_equal(
    query$var_joinids(),
    arrow::concat_arrays(experiment$ms$get("RNA")$var$read()$soma_joinid)
  )

  expect_equal(query$n_obs, n_obs)
  expect_equal(query$n_vars, n_var)

  # X
  expect_error(query$X(), "Must specify a layer name")
  expect_error(query$X(c("a", "b")), "Must specify a single layer name")
  expect_error(query$X("foo"), "The following layer does not exist: foo")

  expect_true(
    query$X("counts")$Equals(
      experiment$ms$get("RNA")$X$get("counts")$read_arrow_table()
    )
  )
})
