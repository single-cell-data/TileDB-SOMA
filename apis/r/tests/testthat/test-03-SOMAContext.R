test_that("SOMAContext with default options", {
  context <- SOMAContext$new()

  # Check config
  config <- context$get_config()
  expect_type(
    sm_ratio <- config["sm.mem.reader.sparse_global_order.ratio_array_data"],
    "character"
  )
  expect_named(sm_ratio)
  expect_equal(unname(sm_ratio), "0.3")

  # Check data protocol for local object
  expect_type(
    data_protocol <- context$get_data_protocol("soma_experiment_local"),
    "character"
  )
  expect_equal(data_protocol, "tiledbv2")

})
