test_that("file path construction handles remote URLs", {
  expect_identical(
    file_path("foo"),
    file.path("foo")
  )
  expect_identical(
    file_path("foo", "bar"),
    file.path("foo", "bar")
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
  expect_equal(uri_scheme("foo/bar"), NULL)
  expect_equal(uri_scheme("/foo/bar"), NULL)
  expect_equal(uri_scheme("file://foo/bar"), "file")
  expect_equal(uri_scheme("file:///foo/bar"), "file")
  expect_equal(uri_scheme("s3://my/bucket"), "s3")
  expect_equal(uri_scheme("tiledb://my/array"), "tiledb")
})

test_that("schemes are removed from uris", {
  expect_equal(uri_scheme_remove("foo/bar"), "foo/bar")
  expect_equal(uri_scheme_remove("/foo/bar"), "/foo/bar")
  expect_equal(uri_scheme_remove("file://foo/bar"), "foo/bar")
  expect_equal(uri_scheme_remove("file:///foo/bar"), "/foo/bar")
  expect_equal(uri_scheme_remove("s3://my/bucket"), "my/bucket")
  expect_equal(uri_scheme_remove("tiledb://my/array"), "my/array")
})

test_that("relative uris are calculated correctly", {
  expect_equal(make_uri_relative("foo/bar", "foo"), "bar")
  expect_equal(make_uri_relative("/foo/bar", "/foo"), "bar")
  expect_equal(make_uri_relative("file://foo/bar", "file://foo"), "bar")

  # Heterogenous schemes
  expect_equal(make_uri_relative("foo/bar", "file://foo"), "bar")
  expect_equal(make_uri_relative("file://foo/bar", "foo"), "bar")

  # Expected errors
  expect_error(
    make_uri_relative("file://bar", "file://foo/bar"),
    "Unable to make relative path between URIs with no common parent"
  )

  expect_error(
    make_uri_relative("s3://foo/bar", "file://foo"),
    "Unable to make relative path between URIs with different schemes"
  )
  expect_error(
    make_uri_relative("s3://foo/bar", "file://foo"),
    "Unable to make relative path between URIs with different schemes"
  )
})
