test_that("SOMAAxisQuery", {
  query <- SOMAAxisQuery$new()
  expect_null(query$value_filter)
  expect_null(query$coords)

  query <- SOMAAxisQuery$new(value_filter = "foo")
  expect_equal(query$value_filter, "foo")
  expect_null(query$coords)

  query <- SOMAAxisQuery$new(coords = list(foo = 1L:2L))
  expect_null(query$value_filter)
  expect_equal(query$coords, list(foo = 1L:2L))

  query <- SOMAAxisQuery$new(value_filter = "foo", coords = list(foo = 1L:2L))
  expect_equal(query$value_filter, "foo")
  expect_equal(query$coords, list(foo = 1L:2L))

  # Bare vector is wrapped in a list
  query <- SOMAAxisQuery$new(coords = 1L:2L)
  expect_equal(query$coords, list(1L:2L))

  # Unnamed list is valid for a single set of coordinates
  query <- SOMAAxisQuery$new(coords = list(1L:2L))
  expect_equal(query$coords, list(1L:2L))


  # Expected failures
  expect_error(
    SOMAAxisQuery$new(value_filter = 1),
    "'value_filter' must be a scalar character"
  )
  expect_error(
    SOMAAxisQuery$new(value_filter = c("foo", "bar")),
    "'value_filter' must be a scalar character"
  )
  expect_error(
    SOMAAxisQuery$new(coords = list(1L:2L, 3L:4L)),
    "'coords' must be a named list to query multiple dimensions"
  )

  expect_error(
    SOMAAxisQuery$new(coords = list(foo = letters)),
    "'coords' must be a list of numeric vectors"
  )
})
