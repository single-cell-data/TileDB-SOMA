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

test_that("returns all coordinates by default", {
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

test_that("querying by dimension coordinates", {
  uri <- withr::local_tempdir("soma-experiment-query-coords")
  n_obs <- 1001L
  n_var <- 99L

  obs_slice <- bit64::as.integer64(seq(3, 72))
  var_slice <- bit64::as.integer64(seq(7, 21))

  experiment <- create_and_populate_experiment(
    uri = uri,
    n_obs = n_obs,
    n_var = n_var,
    X_layer_names = c("counts", "logcounts")
  )

  query <- ExperimentAxisQuery$new(
    experiment = experiment,
    measurement_name = "RNA",
    obs_query = AxisQuery$new(coords = list(soma_joinid = obs_slice)),
    var_query = AxisQuery$new(coords = list(soma_joinid = var_slice))
  )

  expect_true(query$n_obs == diff(range(obs_slice)) + 1)
  expect_true(query$n_vars == diff(range(var_slice)) + 1)

  expect_equal(query$obs()$soma_joinid$as_vector(), as.integer(obs_slice))
  expect_equal(query$var_joinids()$as_vector(), as.integer(var_slice))

  expect_equal(
    query$obs(column_names = "soma_joinid")$soma_joinid$as_vector(),
    as.integer(obs_slice)
  )
  expect_equal(
    query$var(column_names = "soma_joinid")$soma_joinid$as_vector(),
    as.integer(var_slice)
  )

  raw_X <- experiment$ms$get("RNA")$X$get("counts")$read_arrow_table(
    coords = list(obs_slice, var_slice)
  )
  expect_true(query$X("counts")$Equals(raw_X))
})

test_that("querying by value filters", {
  uri <- withr::local_tempdir("soma-experiment-query-value-filters")
  n_obs <- 1001L
  n_var <- 99L

  obs_label_values <- c("1003", "1007", "1038", "1099")
  var_label_values <- c("1018", "1034", "1067")

  experiment <- create_and_populate_experiment(
    uri = uri,
    n_obs = n_obs,
    n_var = n_var,
    X_layer_names = c("counts", "logcounts")
  )

  # TODO: simplify once tiledb-r supports membership expressions
  obs_value_filter <- paste0(
    sprintf("baz == '%s'", obs_label_values),
    collapse = "||"
  )
  var_value_filter <- paste0(
    sprintf("quux == '%s'", var_label_values),
    collapse = "||"
  )

  query <- ExperimentAxisQuery$new(
    experiment = experiment,
    measurement_name = "RNA",
    obs_query = AxisQuery$new(value_filter = obs_value_filter),
    var_query = AxisQuery$new(value_filter = var_value_filter)
  )

  expect_true(query$n_obs == length(obs_label_values))
  expect_true(query$n_vars == length(var_label_values))

  expect_equal(query$obs()$baz$as_vector(), obs_label_values)
  expect_equal(query$var()$quux$as_vector(), var_label_values)
})
