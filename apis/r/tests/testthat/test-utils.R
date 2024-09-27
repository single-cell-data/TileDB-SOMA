test_that("validate read coords", {

  # NULL is a valid value
  expect_equal(
    validate_read_coords(NULL),
    NULL
  )

  # bare vector is converted to list
  expect_equal(
    validate_read_coords(1:10),
    list(1:10)
  )

  expect_equal(
    validate_read_coords(list(1:10)),
    list(1:10)
  )

  # list of multiple coordinate vectors must be named
  expect_error(
    validate_read_coords(list(1:10, 1:10)),
    "'coords' must be a named list to query multiple dimensions"
  )

  expect_equal(
    validate_read_coords(list(int_column = 1:10, float_column = 1:10)),
    list(int_column = 1:10, float_column = 1:10)
  )
})


test_that("validate read coords with dimension names", {

  # assume vector or unnamed list of length 1 corresponds to first dimension
  expect_equal(
    validate_read_coords(1:10, dimnames = "int_column"),
    list(int_column = 1:10)
  )

  # list of named coordinates must match provided dimension names
  expect_error(
    validate_read_coords(list(int_column = 1:10, float_column = 1:10), c("int_column", "string_column")),
    "names of 'coords' must correspond to dimension names"
  )

  expect_error(
    validate_read_coords(list(int_column = 1:10, float_column = 1:10), c("int_column")),
    "names of 'coords' must correspond to dimension names"
  )
})


test_that("validate read coords with dimension names and schema", {
  asch <- create_arrow_schema()

  # if schema is provided, dimnames must be provided too
  expect_error(
    validate_read_coords(1:10, schema = asch),
    "'dimnames' must be provided with a 'schema'"
  )

  # assume vector or unnamed list of length 1 corresponds to first dimension
  expect_equal(
    validate_read_coords(1:10, dimnames = "int_column", schema = asch),
    list(int_column = 1:10)
  )

  # integer coordinates corresponding to int64 dimensions are cast to int64
  expect_equal(
    validate_read_coords(1:10, dimnames = "soma_joinid", schema = asch),
    list(soma_joinid = as.integer64(1:10))
  )

  # casting is selective and only applies to int64 dimensions
  test_coords <- validate_read_coords(
      coords = list(int_column = 1:10, soma_joinid = 1:10),
      dimnames = c("int_column", "soma_joinid"),
      schema = asch
  )

  expect_equal(test_coords$int_column, 1:10)
  expect_equal(test_coords$soma_joinid, as.integer64(1:10))
})
