test_that("SOMAContext with default options", {
  # Create context.
  context <- SOMAContext$new()

  # Check config.
  config <- context$get_config()
  expect_type(
    sm_ratio <- config["sm.mem.reader.sparse_global_order.ratio_array_data"],
    "character"
  )
  expect_named(sm_ratio)
  expect_equal(unname(sm_ratio), "0.3")

  # Check data protocol for local object.
  expect_type(
    data_protocol <- context$get_data_protocol("soma_experiment_local"),
    "character"
  )
  expect_equal(data_protocol, "tiledbv2")

})

test_that("SOMAContext set config", {
  # Create config and context.
  key1 <- "extra.first_key"
  key2 <- "extra.second_key"
  sm_ratio_key <- "sm.mem.reader.sparse_global_order.ratio_array_data"
  expected_value1 <- "9000"
  expected_value2 <- "value2"
  input_config <- c(
    extra.first_key = expected_value1,
    extra.second_key = expected_value2
  )
  context <- SOMAContext$new(input_config)

  # Check config.
  output_config <- context$get_config()
  expect_type(sm_ratio <- output_config[sm_ratio_key], "character")
  expect_equal(unname(sm_ratio), "0.3")
  expect_type(value1 <- output_config[key1], "character")
  expect_equal(unname(value1), expected_value1)
  expect_type(value2 <- output_config[key2], "character")
  expect_equal(unname(value2), expected_value2)
})

test_that("SOMAContext set config with ratio_array_data", {
  # Create config and context.
  key1 <- "extra.first_key"
  key2 <- "extra.second_key"
  sm_ratio_key <- "sm.mem.reader.sparse_global_order.ratio_array_data"
  expected_value1 <- "9000"
  expected_value2 <- "value2"
  expected_sm_ratio <- "0.5"
  input_config <- c(
    sm.mem.reader.sparse_global_order.ratio_array_data = expected_sm_ratio,
    extra.first_key = expected_value1,
    extra.second_key = expected_value2
  )
  context <- SOMAContext$new(input_config)

  # Check config.
  output_config <- context$get_config()
  expect_type(sm_ratio <-output_config[sm_ratio_key], "character")
  expect_equal(unname(sm_ratio), expected_sm_ratio)
  expect_type(value1 <- output_config[key1], "character")
  expect_equal(unname(value1), expected_value1)
  expect_type(value2 <- output_config[key2], "character")
  expect_equal(unname(value2), expected_value2)
})

test_that("SOMAContext from default SOMATileDBContext", {
  # Create context.
  with_mocked_bindings(
    .tiledbsoma_deprecation_version = function() "2.3.0",
    .deprecation_stage = function(when) "deprecate",
    {
      lifecycle::expect_deprecated(tiledbsoma_ctx <- SOMATileDBContext$new())
      lifecycle::expect_deprecated(context <- get_soma_context(NULL, tiledbsoma_ctx, what="get_soma_context()"))
    }
  )

  # Check config.
  config <- context$get_config()
  expect_type(
    sm_ratio <- config["sm.mem.reader.sparse_global_order.ratio_array_data"],
    "character"
  )
  expect_named(sm_ratio)
  expect_equal(unname(sm_ratio), "0.3")

})

test_that("SOMAContext from SOMATileDBContext with config", {
  # Create config and context.
  key1 <- "extra.first_key"
  key2 <- "extra.second_key"
  sm_ratio_key <- "sm.mem.reader.sparse_global_order.ratio_array_data"
  expected_value1 <- "9000"
  expected_value2 <- "value2"
  input_config <- c(
    extra.first_key = expected_value1,
    extra.second_key = expected_value2
  )
  with_mocked_bindings(
    .tiledbsoma_deprecation_version = function() "2.3.0",
    .deprecation_stage = function(when) "deprecate",
    {
      lifecycle::expect_deprecated(tiledbsoma_ctx <- SOMATileDBContext$new(config=input_config))
      lifecycle::expect_deprecated(context <- get_soma_context(NULL, tiledbsoma_ctx, what="get_soma_context()"))
    }
  )

  # Check config.
  output_config <- context$get_config()
  expect_type(sm_ratio <- output_config[sm_ratio_key], "character")
  expect_equal(unname(sm_ratio), "0.3")
  expect_type(value1 <- output_config[key1], "character")
  expect_equal(unname(value1), expected_value1)
  expect_type(value2 <- output_config[key2], "character")
  expect_equal(unname(value2), expected_value2)

})


test_that("SOMAContext from SOMATileDBContext with ratio_array_data", {
  # Create config and context.
  key1 <- "extra.first_key"
  key2 <- "extra.second_key"
  sm_ratio_key <- "sm.mem.reader.sparse_global_order.ratio_array_data"
  expected_value1 <- "9000"
  expected_value2 <- "value2"
  expected_sm_ratio <- "0.5"
  input_config <- c(
    sm.mem.reader.sparse_global_order.ratio_array_data = expected_sm_ratio,
    extra.first_key = expected_value1,
    extra.second_key = expected_value2
  )
  with_mocked_bindings(
    .tiledbsoma_deprecation_version = function() "2.3.0",
    .deprecation_stage = function(when) "deprecate",
    {
      lifecycle::expect_deprecated(tiledbsoma_ctx <- SOMATileDBContext$new(config=input_config))
      lifecycle::expect_deprecated(context <- get_soma_context(NULL, tiledbsoma_ctx, what="get_soma_context(tiledbsoma_ctx)"))
    }
  )

  # Check config.
  output_config <- context$get_config()
  expect_type(sm_ratio <- output_config[sm_ratio_key], "character")
  expect_equal(unname(sm_ratio), "0.5")
  expect_type(value1 <- output_config[key1], "character")
  expect_equal(unname(value1), expected_value1)
  expect_type(value2 <- output_config[key2], "character")
  expect_equal(unname(value2), expected_value2)

})

test_that("Set and get context from environment (deprecated)", {
  # Create config and context.
  key1 <- "extra.first_key"
  key2 <- "extra.second_key"
  sm_ratio_key <- "sm.mem.reader.sparse_global_order.ratio_array_data"
  expected_value1 <- "9000"
  expected_value2 <- "value2"
  expected_sm_ratio <- "0.5"
  input_config <- c(
    sm.mem.reader.sparse_global_order.ratio_array_data = expected_sm_ratio,
    extra.first_key = expected_value1,
    extra.second_key = expected_value2
  )
  with_mocked_bindings(
    .tiledbsoma_deprecation_version = function() "2.3.0",
    .deprecation_stage = function(when) "deprecate",
    {
       lifecycle::expect_deprecated(tiledbsoma::soma_context(config=input_config))
    }
  )
  context <- get_soma_context(NULL, NULL, what="get_soma_context(tiledbsoma_ctx)")

  # Check config.
  output_config <- context$get_config()
  expect_type(sm_ratio <- output_config[sm_ratio_key], "character")
  expect_equal(unname(sm_ratio), "0.5")
  expect_type(value1 <- output_config[key1], "character")
  expect_equal(unname(value1), expected_value1)
  expect_type(value2 <- output_config[key2], "character")
  expect_equal(unname(value2), expected_value2)
})
