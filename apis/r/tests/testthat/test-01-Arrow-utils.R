test_that("Validating arrow data type compatibility", {
  expect_false(check_arrow_data_types(arrow::int32(), arrow::float32()))
  expect_true(check_arrow_data_types(arrow::int32(), arrow::int32()))
  # strings and large strings are compatible
  expect_true(check_arrow_data_types(arrow::string(), arrow::large_utf8()))

  expect_error(
    check_arrow_data_types("not an arrow data type", arrow::int32())
  )
})

test_that("Validating arrow schema data type compatibility", {
  from <- arrow::schema(int_column = arrow::int32())
  to <- arrow::schema(int_column = arrow::int32())
  expect_true(check_arrow_schema_data_types(from, to))

  # Add incompatible fields
  from$float_column <- arrow::int16()
  to$float_column <- arrow::float16()
  expect_error(
    check_arrow_schema_data_types(from, to),
    "Schemas are incompatible"
  )

  # Schemas with different fields
  from$string_column <- arrow::string()
  expect_error(
    check_arrow_schema_data_types(from, to),
    "'from' and 'to' must have the same number of fields"
  )

  to$fizz <- arrow::string()
  expect_error(
    check_arrow_schema_data_types(from, to),
    "'from' and 'to' must have the same field names"
  )
})
