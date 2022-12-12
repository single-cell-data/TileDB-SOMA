
test_that("SOMACollection basics", {
  uri <- file.path(withr::local_tempdir(), "new-collection")
  collection <- SOMACollection$new(uri)
  expect_equal(collection$uri, uri)

  # Should not exist on disk until created
  expect_false(dir.exists(uri))
  expect_false(collection$exists())

  # Create the collection on disk
  collection$create()

  expect_true(dir.exists(uri))
  expect_match(tiledb::tiledb_object_type(uri), "GROUP")
  expect_true(collection$soma_type == "SOMACollection")
  expect_true(collection$exists())
  expect_equal(collection$length(), 0)

  # Add an element to the collection
  dataframe <- create_and_populate_soma_dataframe(file.path(uri, "sdf"))
  collection$set(dataframe, name = "sdf")

  # Read back the collection
  readback_collection <- SOMACollection$new(uri)
  expect_equal(readback_collection$length(), 1)

  readback_dataframe <- readback_collection$get("sdf")
  expect_is(readback_dataframe, "SOMADataFrame")
})
