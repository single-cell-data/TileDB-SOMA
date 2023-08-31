test_that("TileDB classes can be converted to Arrow equivalents", {

  # Dimension to Arrow field
  dim0 <- tiledb::tiledb_dim(
    name = "dim0",
    domain = NULL,
    tile = NULL,
    type = "ASCII"
  )

  dim1 <- tiledb::tiledb_dim(
    name = "dim1",
    domain = bit64::as.integer64(c(0, 100)),
    tile = bit64::as.integer64(10),
    type = "INT64"
  )

  # Error if not a dimension
  expect_error(arrow_field_from_tiledb_attr(dim0))

  # String dimension
  dim0_field <- arrow_field_from_tiledb_dim(dim0)
  expect_true(is_arrow_field(dim0_field))
  expect_equal(dim0_field$name, tiledb::name(dim0))
  expect_equal(
    tiledb_type_from_arrow_type(dim0_field$type),
    tiledb::datatype(dim0)
  )

  # Integer dimension
  dim1_field <- arrow_field_from_tiledb_dim(dim1)
  expect_true(is_arrow_field(dim1_field))
  expect_equal(dim1_field$name, tiledb::name(dim1))
  expect_equal(
    tiledb_type_from_arrow_type(dim1_field$type),
    tiledb::datatype(dim1)
  )

  # Attribute to Arrow field
  attr0 <- tiledb::tiledb_attr(
    name = "attr0",
    type = "ASCII"
  )

  attr1 <- tiledb::tiledb_attr(
    name = "attr1",
    type = "INT64"
  )

  # Error if not an attribute
  expect_error(arrow_field_from_tiledb_attr(dim0))

  # String attribute
  attr0_field <- arrow_field_from_tiledb_attr(attr0)
  expect_true(is_arrow_field(attr0_field))
  expect_equal(attr0_field$name, tiledb::name(attr0))
  expect_equal(
    tiledb_type_from_arrow_type(attr0_field$type),
    tiledb::datatype(attr0)
  )

  # Integer attribute
  attr1_field <- arrow_field_from_tiledb_attr(attr1)
  expect_true(is_arrow_field(attr1_field))
  expect_equal(attr1_field$name, tiledb::name(attr1))
  expect_equal(
    tiledb_type_from_arrow_type(attr1_field$type),
    tiledb::datatype(attr1)
  )

  # TileDB schema to Arrow schema
  tdb_schema <- tiledb::tiledb_array_schema(
    domain = tiledb::tiledb_domain(c(dim0, dim1)),
    attrs = c(attr0, attr1),
    sparse = TRUE
  )

  arrow_schema <- arrow_schema_from_tiledb_schema(tdb_schema)
  expect_true(is_arrow_schema(arrow_schema))
  expect_equal(length(arrow_schema$fields), 4)
  expect_equal(names(arrow_schema), c("dim0", "dim1", "attr0", "attr1"))
})

# #1651 CI debug start
#test_that("Validating arrow data type compatibility", {
#  expect_false(check_arrow_data_types(arrow::int32(), arrow::float32()))
#  expect_true(check_arrow_data_types(arrow::int32(), arrow::int32()))
#  # strings and large strings are compatible
#  expect_true(check_arrow_data_types(arrow::string(), arrow::large_utf8()))
#
#  expect_error(
#    check_arrow_data_types("not an arrow data type", arrow::int32())
#  )
#})
#
#test_that("Validating arrow schema data type compatibility", {
#  from <- arrow::schema(foo = arrow::int32())
#  to <- arrow::schema(foo = arrow::int32())
#  expect_true(check_arrow_schema_data_types(from, to))
#
#  # Add incompatible fields
#  from$bar <- arrow::int16()
#  to$bar <- arrow::float16()
#  expect_error(
#    check_arrow_schema_data_types(from, to),
#    "Schemas are incompatible"
#  )
#
#  # Schemas with different fields
#  from$baz <- arrow::string()
#  expect_error(
#    check_arrow_schema_data_types(from, to),
#    "'from' and 'to' must have the same number of fields"
#  )
#
#  to$fizz <- arrow::string()
#  expect_error(
#    check_arrow_schema_data_types(from, to),
#    "'from' and 'to' must have the same field names"
#  )
#})
