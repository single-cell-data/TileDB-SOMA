test_that("TileDBArray helper functions", {
  uri <- withr::local_tempdir(pattern = "test-array")

  tdba <- TileDBArray$new(uri = uri, internal_use_only = "allowed_use")

  # Check errors on non-existent array
  expect_error(tdba$object, "Array does not exist.")
  ## inactive under new implementation  expect_error(tdba$set_metadata(list(foo = "bar")), "Item must be open for write.")

  # Create an array
  index_cols <- c("Dept", "Gender")
  df <- as.data.frame(UCBAdmissions)
  tiledb::fromDataFrame(df, uri, col_index = index_cols)

  expect_identical(tdba$uri, uri)
  expect_is(tdba$tiledb_array(), "tiledb_array")
  expect_is(tdba$object, "tiledb_array")
  expect_identical(tdba$dimnames(), index_cols)
  expect_identical(tdba$index_column_names(), index_cols)

  attr_cols <- setdiff(colnames(df), index_cols)
  expect_identical(tdba$attrnames(), attr_cols)

  # both dimensions and attributes
  expect_setequal(tdba$colnames(), colnames(df))

  # metadata
  md <- list(baz = "qux", foo = "bar")
  tdba$open(mode = "WRITE", internal_use_only = "allowed_use")
  tdba$set_metadata(md)
  tdba$close()

  tdba$open(mode = "READ", internal_use_only = "allowed_use")
  expect_equal(tdba$get_metadata(key = "foo"), "bar")
  expect_equal(tdba$get_metadata(key = "baz"), "qux")
  expect_equal(length(tdba$get_metadata()), 2)
  tdba$close()

  # The SOMA spec requires the ability to read back metadata even when the
  # array is opened for write.
  tdba$open(mode = "WRITE", internal_use_only = "allowed_use")
  expect_equal(tdba$get_metadata(key = "foo"), "bar")
  expect_equal(tdba$get_metadata(key = "baz"), "qux")
  expect_equal(length(tdba$get_metadata()), 2)
  tdba$close()

  ## shape
  tdba$open(mode = "READ", internal_use_only = "allowed_use")
  expect_equal(tdba$ndim(), 2)
  tdba$close()
})
