test_that("local object uris using file:// scheme are detected", {
  uri <- "/data/object_name"
  tdb_uri <- TileDBURI$new(uri)

  expect_match(tdb_uri$uri, uri)
  expect_false(tdb_uri$is_remote_uri())
  expect_false(tdb_uri$is_tiledb_cloud_uri())
  expect_false(tdb_uri$is_tiledb_cloud_creation_uri())
})

test_that("local object uris using file:// scheme are detected", {
  uri <- "file://data/object_name"
  tdb_uri <- TileDBURI$new(uri)

  expect_match(tdb_uri$uri, uri)
  expect_false(tdb_uri$is_remote_uri())
  expect_false(tdb_uri$is_tiledb_cloud_uri())
  expect_false(tdb_uri$is_tiledb_cloud_creation_uri())
})

test_that("s3 object uris are handled", {
  uri <- "s3://bucket/prefix"
  tdb_uri <- TileDBURI$new(uri)

  expect_false(tdb_uri$is_tiledb_cloud_uri())
  expect_null(tdb_uri$tiledb_cloud_uri)

  expect_false(tdb_uri$is_tiledb_cloud_creation_uri())
  expect_equal(tdb_uri$object_uri, uri)
})

test_that("tiledb cloud uris are handled", {
  uri <- "tiledb://namespace/object_name"
  tdb_uri <- TileDBURI$new(uri)

  expect_true(tdb_uri$is_tiledb_cloud_uri())
  expect_equal(tdb_uri$tiledb_cloud_uri, uri)

  expect_false(tdb_uri$is_tiledb_cloud_creation_uri())
  expect_equal(tdb_uri$object_uri, uri)
})

test_that("tiledb cloud creation uris are detected", {
  uri <- "tiledb://namespace/s3://bucket/prefix/object_name"
  tdb_uri <- TileDBURI$new(uri)
  expect_equal(tdb_uri$uri, uri)

  expect_true(tdb_uri$is_tiledb_cloud_uri())
  expect_equal(tdb_uri$tiledb_cloud_uri, "tiledb://namespace/object_name")

  expect_true(tdb_uri$is_tiledb_cloud_creation_uri())
  expect_equal(tdb_uri$object_uri, "s3://bucket/prefix/object_name")
})
