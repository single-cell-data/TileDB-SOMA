test_that("TileDBArray helper functions", {
  uri <- withr::local_tempdir(pattern = "test-array")

  tdba <- TileDBArray$new(uri = uri, internal_use_only = "allowed_use")

  # Check errors on non-existent array
  expect_error(tdba$object, "Array does not exist.")
  expect_error(tdba$set_metadata(list(foo = "bar")), "Item must be open for write.")

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

  tdba$open(mode = "READ", internal_use_only = "allowed_use")

  # dimension slicing
  expect_error(
    tdba$set_query(dims = "foo"),
    "'dims' must be a named list of character vectors"
  )
  expect_error(
    tdba$set_query(dims = "foo"),
    "'dims' must be a named list of character vectors"
  )
  expect_error(
    tdba$set_query(dims = list(a = 1L)),
    "'dims' must be a named list of character vectors"
  )
  expect_error(
    tdba$set_query(dims = list(foo = "bar")),
    "The following dimension does not exist: foo"
  )

  expect_silent(
    tdba$set_query(dims = list(Dept = c("A", "B")))
  )

  # verify selected ranges were set
  expect_equal(
    tiledb::selected_ranges(tdba$object),
    list(Dept = cbind(c("A", "B"), c("A", "B")))
  )

  # query result includes only selected ranges
  expect_equal(
    unique(tdba$object[]$Dept),
    c("A", "B")
  )
  tdba$close()

  # set attribute filter
  tdba$open(mode = "READ", internal_use_only = "allowed_use")
  tdba$set_query(attr_filter = Admit == "Admitted")
  expect_true(all(tdba$object[]$Admit == "Admitted"))

  # update attribute filter
  tdba$set_query(attr_filter = Admit != "Admitted")
  expect_true(all(tdba$object[]$Admit == "Rejected"))

  # reset attribute filter
  tdba$reset_query()
  expect_length(tdba$object[]$Admit, nrow(df))

  ## shape
  expect_equal(tdba$ndim(), 2)

  tdba$close()
})
