test_that("file path construction handles remote URLs", {
  expect_identical(
    file_path("int_column"),
    file.path("int_column")
  )
  expect_identical(
    file_path("int_column", "bar"),
    file.path("int_column", "bar")
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
  expect_equal(uri_scheme("int_column/bar"), NULL)
  expect_equal(uri_scheme("/int_column/bar"), NULL)
  expect_equal(uri_scheme("file://int_column/bar"), "file")
  expect_equal(uri_scheme("file:///int_column/bar"), "file")
  expect_equal(uri_scheme("s3://my/bucket"), "s3")
  expect_equal(uri_scheme("tiledb://my/array"), "tiledb")
})

test_that("schemes are removed from uris", {
  expect_equal(uri_scheme_remove("int_column/bar"), "int_column/bar")
  expect_equal(uri_scheme_remove("/int_column/bar"), "/int_column/bar")
  expect_equal(uri_scheme_remove("file://int_column/bar"), "int_column/bar")
  expect_equal(uri_scheme_remove("file:///int_column/bar"), "/int_column/bar")
  expect_equal(uri_scheme_remove("s3://my/bucket"), "my/bucket")
  expect_equal(uri_scheme_remove("tiledb://my/array"), "my/array")
})

test_that("relative uris are calculated correctly", {
  expect_equal(make_uri_relative("int_column/bar", "int_column"), "bar")
  expect_equal(make_uri_relative("/int_column/bar", "/int_column"), "bar")
  expect_equal(make_uri_relative("file://int_column/bar", "file://int_column"), "bar")

  # Heterogenous schemes
  expect_equal(make_uri_relative("int_column/bar", "file://int_column"), "bar")
  expect_equal(make_uri_relative("file://int_column/bar", "int_column"), "bar")

  # Expected errors
  expect_error(
    make_uri_relative("file://bar", "file://int_column/bar"),
    "Unable to make relative path between URIs with no common parent"
  )

  expect_error(
    make_uri_relative("s3://int_column/bar", "file://int_column"),
    "Unable to make relative path between URIs with different schemes"
  )
  expect_error(
    make_uri_relative("s3://int_column/bar", "file://int_column"),
    "Unable to make relative path between URIs with different schemes"
  )
})
