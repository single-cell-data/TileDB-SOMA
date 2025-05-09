test_that("returns all coordinates by default", {
  skip_if(!extended_tests())
  uri <- tempfile(pattern = "soma-experiment-query-all")
  n_obs <- 20L
  n_var <- 10L

  experiment <- create_and_populate_experiment(
    uri = uri,
    n_obs = n_obs,
    n_var = n_var,
    X_layer_names = c("counts", "logcounts"),
    mode = "READ"
  )
  on.exit(experiment$close())

  query <- SOMAExperimentAxisQuery$new(
    experiment = experiment,
    measurement_name = "RNA"
  )

  # Query object can be created from the SOMAExperiment axis_query() method
  expect_equal(experiment$axis_query(measurement_name = "RNA"), query)

  # obs/var tables
  expect_true(is(query$obs(), "TableReadIter"))
  expect_true(is(query$var(), "TableReadIter"))
  expect_true(query$obs()$concat()$Equals(experiment$obs$read()$concat()))
  expect_true(query$var()$concat()$Equals(experiment$ms$get("RNA")$var$read()$concat()))

  # obs/var joinids
  expect_equal(
    query$obs_joinids(),
    arrow::concat_arrays(experiment$obs$read()$concat()$soma_joinid)
  )
  expect_equal(
    query$var_joinids(),
    arrow::concat_arrays(experiment$ms$get("RNA")$var$read()$concat()$soma_joinid)
  )

  expect_equal(query$n_obs, n_obs)
  expect_equal(query$n_vars, n_var)

  # X
  expect_error(query$X(), "Must specify an X layer name")
  expect_error(query$X(c("a", "b")), "Must specify a single X layer name")
  expect_error(query$X("int_column"), "The following layer does not exist: int_column")

  expect_true(is(query$X("counts"), "SOMASparseNDArrayRead"))
  expect_true(
    query$X("counts")$tables()$concat()$Equals(
      experiment$ms$get("RNA")$X$get("counts")$read()$tables()$concat()
    )
  )

  experiment$close()
})

test_that("querying by dimension coordinates", {
  skip_if(!extended_tests())
  uri <- tempfile(pattern = "soma-experiment-query-coords")
  n_obs <- 1001L
  n_var <- 99L

  obs_slice <- bit64::as.integer64(seq(3, 72))
  var_slice <- bit64::as.integer64(seq(7, 21))

  experiment <- create_and_populate_experiment(
    uri = uri,
    n_obs = n_obs,
    n_var = n_var,
    X_layer_names = c("counts", "logcounts"),
    mode = "READ"
  )
  on.exit(experiment$close())

  query <- SOMAExperimentAxisQuery$new(
    experiment = experiment,
    measurement_name = "RNA",
    obs_query = SOMAAxisQuery$new(coords = list(soma_joinid = obs_slice)),
    var_query = SOMAAxisQuery$new(coords = list(soma_joinid = var_slice))
  )

  expect_true(query$n_obs == diff(range(obs_slice)) + 1)
  expect_true(query$n_vars == diff(range(var_slice)) + 1)

  expect_equal(query$obs()$concat()$soma_joinid$as_vector(), as.integer(obs_slice))
  expect_equal(query$var_joinids()$as_vector(), as.integer(var_slice))

  expect_equal(
    query$obs(column_names = "soma_joinid")$concat()$soma_joinid$as_vector(),
    as.integer(obs_slice)
  )
  expect_equal(
    query$var(column_names = "soma_joinid")$concat()$soma_joinid$as_vector(),
    as.integer(var_slice)
  )

  raw_X <- experiment$ms$get("RNA")$X$get("counts")$read(
    coords = list(obs_slice, var_slice)
  )$tables()$concat()
  expect_true(query$X("counts")$tables()$concat()$Equals(raw_X))

  experiment$close()
})

test_that("querying by value filters", {
  skip_if(!extended_tests())
  uri <- tempfile(pattern = "soma-experiment-query-value-filters")
  n_obs <- 1001L
  n_var <- 99L

  obs_label_values <- c("1003", "1007", "1038", "1099")
  var_label_values <- c("1018", "1034", "1067")

  experiment <- create_and_populate_experiment(
    uri = uri,
    n_obs = n_obs,
    n_var = n_var,
    X_layer_names = c("counts", "logcounts"),
    mode = "READ"
  )
  on.exit(experiment$close())

  # TODO: simplify once tiledb-r supports membership expressions
  obs_value_filter <- paste0(
    sprintf("string_column == '%s'", obs_label_values),
    collapse = "||"
  )
  var_value_filter <- paste0(
    sprintf("quux == '%s'", var_label_values),
    collapse = "||"
  )

  query <- SOMAExperimentAxisQuery$new(
    experiment = experiment,
    measurement_name = "RNA",
    obs_query = SOMAAxisQuery$new(value_filter = obs_value_filter),
    var_query = SOMAAxisQuery$new(value_filter = var_value_filter)
  )

  expect_true(query$n_obs == length(obs_label_values))
  expect_true(query$n_vars == length(var_label_values))

  expect_equal(query$obs()$concat()$string_column$as_vector(), obs_label_values)
  expect_equal(query$var()$concat()$quux$as_vector(), var_label_values)

  experiment$close()
})

test_that("query by value filters with enums", {
  skip_if(!extended_tests())
  uri <- tempfile(pattern = "soma-experiment-query-enum-filters")
  n_obs <- 1001L
  n_var <- 99L

  obs_label_values <- c("1003", "1007", "1038", "1099")
  var_label_values <- c("1018", "1034", "1067")

  experiment <- create_and_populate_experiment(
    uri = uri,
    n_obs = n_obs,
    n_var = n_var,
    X_layer_names = c("counts", "logcounts"),
    mode = "READ"
  )
  on.exit(experiment$close())

  obs_tbl <- experiment$obs$read()$concat()
  obs_tbl$enum <- factor(
    sample(c("red", "blue", "green"), size = n_obs, replace = TRUE),
    levels = c("red", "blue", "green")
  )

  experiment$close()
  experiment <- SOMAExperimentOpen(experiment$uri, mode = "WRITE")
  experiment$update_obs(obs_tbl)

  experiment$close()
  experiment <- SOMAExperimentOpen(experiment$uri, mode = "READ")

  # Test enum query with present level
  query <- SOMAExperimentAxisQuery$new(
    experiment = experiment,
    measurement_name = "RNA",
    obs_query = SOMAAxisQuery$new(value_filter = "enum == 'green'")
  )
  expect_equal(query$n_obs, sum(as.vector(obs_tbl$enum$as_vector()) == "green"))
  obs_df <- as.data.frame(query$obs()$concat())
  expect_s3_class(obs_df$enum, "factor")
  expect_identical(levels(obs_df$enum), c("red", "blue", "green"))
  expect_identical(unique(as.vector(obs_df$enum)), "green")

  # Test enum query with present and missing level
  core <- list(
    tiledbsoma = numeric_version(tiledbsoma:::libtiledbsoma_version(TRUE)),
    tiledb.r = numeric_version(paste(get_tiledb_version(), collapse = "."))
  )
  skip_if(
    any(vapply(core, \(x) x < "2.21", FUN.VALUE = logical(1L))),
    message = "Handling of missing enum levels is implemented in Core 2.21 and higher"
  )

  query <- SOMAExperimentAxisQuery$new(
    experiment = experiment,
    measurement_name = "RNA",
    obs_query = SOMAAxisQuery$new(value_filter = "enum == 'green' || enum == 'purple'")
  )
  expect_equal(query$n_obs, sum(as.vector(obs_tbl$enum$as_vector()) == "green"))
  obs_df <- as.data.frame(query$obs()$concat())
  expect_s3_class(obs_df$enum, "factor")
  expect_identical(levels(obs_df$enum), c("red", "blue", "green"))
  expect_identical(unique(as.vector(obs_df$enum)), "green")

  # Test enum query for everything except missing level
  query <- SOMAExperimentAxisQuery$new(
    experiment = experiment,
    measurement_name = "RNA",
    obs_query = SOMAAxisQuery$new(value_filter = "enum != 'purple'")
  )
  expect_equal(query$n_obs, n_obs)
  obs_df <- as.data.frame(query$obs()$concat())
  expect_equal(nrow(obs_df), n_obs)

  # Test enum query with missing level
  query <- SOMAExperimentAxisQuery$new(
    experiment = experiment,
    measurement_name = "RNA",
    obs_query = SOMAAxisQuery$new(value_filter = "enum == 'purple'")
  )
  expect_equal(query$n_obs, 0L)
})

test_that("querying by both coordinates and value filters", {
  skip_if(!extended_tests())
  uri <- tempfile(pattern = "soma-experiment-query-coords-and-value-filters")

  n_obs <- 1001L
  n_var <- 99L

  obs_slice <- bit64::as.integer64(seq(3, 72))
  var_slice <- bit64::as.integer64(seq(7, 21))

  obs_label_values <- c("1003", "1007", "1038", "1099")
  var_label_values <- c("1018", "1034", "1067")

  # TODO: simplify once tiledb-r supports membership expressions
  obs_value_filter <- paste0(
    sprintf("string_column == '%s'", obs_label_values),
    collapse = "||"
  )
  var_value_filter <- paste0(
    sprintf("quux == '%s'", var_label_values),
    collapse = "||"
  )

  experiment <- create_and_populate_experiment(
    uri = uri,
    n_obs = n_obs,
    n_var = n_var,
    X_layer_names = c("counts", "logcounts"),
    mode = "READ"
  )
  on.exit(experiment$close())

  # obs slice / var value filter
  query <- SOMAExperimentAxisQuery$new(
    experiment = experiment,
    measurement_name = "RNA",
    obs_query = SOMAAxisQuery$new(
      coords = list(soma_joinid = obs_slice),
    ),
    var_query = SOMAAxisQuery$new(
      value_filter = var_value_filter
    )
  )

  expect_true(query$n_obs == diff(range(obs_slice)) + 1)
  expect_true(query$n_vars == length(var_label_values))

  #  obs value filter / var slice
  query <- SOMAExperimentAxisQuery$new(
    experiment = experiment,
    measurement_name = "RNA",
    obs_query = SOMAAxisQuery$new(
      value_filter = obs_value_filter
    ),
    var_query = SOMAAxisQuery$new(
      coords = list(soma_joinid = var_slice)
    )
  )

  expect_true(query$n_obs == length(obs_label_values))
  expect_true(query$n_vars == diff(range(var_slice)) + 1)

  # obs/var slice and value filter
  query <- SOMAExperimentAxisQuery$new(
    experiment = experiment,
    measurement_name = "RNA",
    obs_query = SOMAAxisQuery$new(
      coords = list(soma_joinid = obs_slice),
      value_filter = obs_value_filter
    ),
    var_query = SOMAAxisQuery$new(
      coords = list(soma_joinid = var_slice),
      value_filter = var_value_filter
    )
  )

  # Determine expected results
  obs_df <- experiment$obs$read()$concat()$to_data_frame()
  obs_hits <- obs_df$soma_joinid %in% as.integer(obs_slice) &
    obs_df$string_column %in% obs_label_values

  var_df <- experiment$ms$get("RNA")$var$read()$concat()$to_data_frame()
  var_hits <- var_df$soma_joinid %in% as.integer(var_slice) &
    var_df$quux %in% var_label_values

  expect_equivalent(query$obs()$concat()$to_data_frame(), obs_df[obs_hits, ])
  expect_equivalent(query$var()$concat()$to_data_frame(), var_df[var_hits, ])

  experiment$close()
})

test_that("queries with empty results", {
  skip_if(!extended_tests())
  uri <- tempfile(pattern = "soma-experiment-query-empty-results")
  n_obs <- 1001L
  n_var <- 99L

  experiment <- create_and_populate_experiment(
    uri = uri,
    n_obs = n_obs,
    n_var = n_var,
    X_layer_names = c("counts", "logcounts"),
    mode = "READ"
  )
  on.exit(experiment$close())

  # obs/var slice and value filter
  query <- SOMAExperimentAxisQuery$new(
    experiment = experiment,
    measurement_name = "RNA",
    obs_query = SOMAAxisQuery$new(
      value_filter = "string_column == 'does-not-exist'"
    ),
    var_query = SOMAAxisQuery$new(
      value_filter = "quux == 'does-not-exist'"
    )
  )
  expect_equal(query$obs()$concat()$num_rows, 0)
  expect_equal(query$var()$concat()$num_rows, 0)
})

test_that("retrieving query results in supported formats", {
  skip_if(!extended_tests())
  uri <- tempfile(pattern = "soma-experiment-query-results-formats1")
  n_obs <- 1001L
  n_var <- 99L

  experiment <- create_and_populate_experiment(
    uri = uri,
    n_obs = n_obs,
    n_var = n_var,
    X_layer_names = c("counts", "logcounts"),
    mode = "READ"
  )
  on.exit(experiment$close())

  query <- SOMAExperimentAxisQuery$new(
    experiment = experiment,
    measurement_name = "RNA"
  )

  # Check SOMAAxisQueryResult class
  res <- query$read()
  expect_true(inherits(res, "SOMAAxisQueryResult"))
  expect_true(is_arrow_table(res$obs))
  expect_true(is_arrow_table(res$var))
  expect_true(is.list(res$X_layers))
  expect_true(is_arrow_table(res$X_layers[[1]]))
  expect_true(is_arrow_table(res$X_layers[[2]]))

  experiment$close()
})

test_that("query result value indexer", {
  skip_if(!extended_tests())
  uri <- tempfile(pattern = "soma-experiment-query-results-indexer")
  n_obs <- 1001L
  n_var <- 99L

  obs_slice <- bit64::as.integer64(seq(1, 10))
  var_slice <- bit64::as.integer64(seq(1, 10))

  experiment <- create_and_populate_experiment(
    uri = uri,
    n_obs = n_obs,
    n_var = n_var,
    X_layer_names = c("counts", "logcounts"),
    mode = "READ"
  )
  on.exit(experiment$close())

  query <- SOMAExperimentAxisQuery$new(
    experiment = experiment,
    measurement_name = "RNA",
    obs_query = SOMAAxisQuery$new(
      coords = list(soma_joinid = obs_slice)
    ),
    var_query = SOMAAxisQuery$new(
      coords = list(soma_joinid = var_slice)
    )
  )

  indexer <- query$indexer

  # coords inside query result are indexed
  expect_equal(
    indexer$by_obs(c(1, 4, 2)),
    arrow::Array$create(c(0, 3, 1), type = arrow::int32())
  )

  expect_equal(
    indexer$by_var(c(10, 1)),
    arrow::Array$create(c(9, 0), type = arrow::int32())
  )

  # coords outside query result return null
  expect_equal(
    indexer$by_obs(c(1, 4, 2, 1000)),
    arrow::Array$create(c(0, 3, 1, NA), type = arrow::int32())
  )

  expect_equal(
    indexer$by_var(c(10, 1, 1000)),
    arrow::Array$create(c(9, 0, NA), type = arrow::int32())
  )

  # coord validation
  expect_error(
    indexer$by_obs(),
    "argument \"coords\" is missing, with no default"
  )

  expect_error(
    indexer$by_obs(c(1, 4, 2, 1000, "int_column")),
    "'coords' must be a numeric vector or arrow Array"
  )

  experiment$close()
})

test_that("query result value indexer upcast", {
  skip_if(!extended_tests())
  uri <- tempfile(pattern = "soma-experiment-query-results-indexer-upcast")
  n_obs <- 1001L
  n_var <- 99L

  obs_slice <- seq(1, 10)
  var_slice <- seq(1, 10)

  experiment <- create_and_populate_experiment(
    uri = uri,
    n_obs = n_obs,
    n_var = n_var,
    X_layer_names = c("counts", "logcounts"),
    mode = "READ"
  )
  on.exit(experiment$close())

  query <- SOMAExperimentAxisQuery$new(
    experiment = experiment,
    measurement_name = "RNA",
    obs_query = SOMAAxisQuery$new(
      coords = list(soma_joinid = obs_slice)
    ),
    var_query = SOMAAxisQuery$new(
      coords = list(soma_joinid = var_slice)
    )
  )

  indexer <- query$indexer

  # coords inside query result are indexed
  expect_equal(
    indexer$by_obs(c(1, 4, 2)),
    arrow::Array$create(c(0, 3, 1), type = arrow::int32())
  )

  expect_equal(
    indexer$by_var(c(10, 1)),
    arrow::Array$create(c(9, 0), type = arrow::int32())
  )

  # coords outside query result return null
  expect_equal(
    indexer$by_obs(c(1, 4, 2, 1000)),
    arrow::Array$create(c(0, 3, 1, NA), type = arrow::int32())
  )

  expect_equal(
    indexer$by_var(c(10, 1, 1000)),
    arrow::Array$create(c(9, 0, NA), type = arrow::int32())
  )

  # coord validation
  expect_error(
    indexer$by_obs(),
    "argument \"coords\" is missing, with no default"
  )

  expect_error(
    indexer$by_obs(c(1, 4, 2, 1000, "int_column")),
    "'coords' must be a numeric vector or arrow Array"
  )

  experiment$close()
})
