test_that("file path construction handles remote URLs", {
  expect_identical(
    file_path("int_column"),
    file.path("int_column")
  )
  expect_identical(
    file_path("int_column", "float_column"),
    file.path("int_column", "float_column")
  )
  expect_identical(
    file_path("s3://my", "bucket", fsep = "\\"),
    "s3://my/bucket"
  )
  expect_identical(
    file_path("tiledb://my", "array", fsep = "\\"),
    "tiledb://my/array"
  )
})

test_that("uri schemes are retrieved", {
  expect_equal(uri_scheme("int_column/float_column"), NULL)
  expect_equal(uri_scheme("/int_column/float_column"), NULL)
  expect_equal(uri_scheme("file://int_column/float_column"), "file")
  expect_equal(uri_scheme("file:///int_column/float_column"), "file")
  expect_equal(uri_scheme("s3://my/bucket"), "s3")
  expect_equal(uri_scheme("tiledb://my/array"), "tiledb")
})

test_that("schemes are removed from uris", {
  expect_equal(
    uri_scheme_remove("int_column/float_column"),
    "int_column/float_column"
  )
  expect_equal(
    uri_scheme_remove("/int_column/float_column"),
    "/int_column/float_column"
  )
  expect_equal(
    uri_scheme_remove("file://int_column/float_column"),
    "int_column/float_column"
  )
  expect_equal(
    uri_scheme_remove("file:///int_column/float_column"),
    "/int_column/float_column"
  )
  expect_equal(uri_scheme_remove("s3://my/bucket"), "my/bucket")
  expect_equal(uri_scheme_remove("tiledb://my/array"), "my/array")
})

test_that("relative uris are calculated correctly", {
  expect_equal(
    make_uri_relative("int_column/float_column", "int_column"),
    "float_column"
  )
  expect_equal(
    make_uri_relative("/int_column/float_column", "/int_column"),
    "float_column"
  )
  expect_equal(
    make_uri_relative("file://int_column/float_column", "file://int_column"),
    "float_column"
  )

  # Heterogenous schemes
  expect_equal(
    make_uri_relative("int_column/float_column", "file://int_column"),
    "float_column"
  )
  expect_equal(
    make_uri_relative("file://int_column/float_column", "int_column"),
    "float_column"
  )

  # Expected errors
  expect_error(
    make_uri_relative("file://float_column", "file://int_column/float_column"),
    "Unable to make relative path between URIs with no common parent"
  )

  expect_error(
    make_uri_relative("s3://int_column/float_column", "file://int_column"),
    "Unable to make relative path between URIs with different schemes"
  )
  expect_error(
    make_uri_relative("s3://int_column/float_column", "file://int_column"),
    "Unable to make relative path between URIs with different schemes"
  )
})
