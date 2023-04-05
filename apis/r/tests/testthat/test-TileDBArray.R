test_that("TileDBArray helper functions", {
  uri <- withr::local_tempdir(pattern = "test-array")
  tdb <- TileDBArray$new(uri = uri, internal_use_only = "allowed_use")

  # Check errors on non-existent array
  expect_error(tdb$object, "Array does not exist.")
  expect_error(tdb$set_metadata(list(foo = "bar")), "Array does not exist.")

  # create an array
  index_cols <- c("Dept", "Gender")
  df <- as.data.frame(UCBAdmissions)
  tiledb::fromDataFrame(df, uri, col_index = index_cols)

  expect_identical(tdb$uri, uri)
  expect_is(tdb$tiledb_array(), "tiledb_array")
  expect_is(tdb$object, "tiledb_array")
  expect_identical(tdb$dimnames(), index_cols)
  expect_identical(tdb$index_column_names(), index_cols)

  attr_cols <- setdiff(colnames(df), index_cols)
  expect_identical(tdb$attrnames(), attr_cols)

  # both dimensions and attributes
  expect_setequal(tdb$colnames(), colnames(df))

  # metadata
  md <- list(baz = "qux", foo = "bar")
  tdb$set_metadata(md)
  expect_equal(tdb$get_metadata(key = "foo"), "bar")
  expect_equal(tdb$get_metadata(prefix = "foo"), md["foo"])
  expect_equal(tdb$get_metadata(), md)

  # dimension slicing
  tdb <- TileDBArray$new(uri = uri, internal_use_only = "allowed_use")
  expect_error(
    tdb$set_query(dims = "foo"),
    "'dims' must be a named list of character vectors"
  )
  expect_error(
    tdb$set_query(dims = "foo"),
    "'dims' must be a named list of character vectors"
  )
  expect_error(
    tdb$set_query(dims = list(a = 1L)),
    "'dims' must be a named list of character vectors"
  )
  expect_error(
    tdb$set_query(dims = list(foo = "bar")),
    "The following dimension does not exist: foo"
  )

  expect_silent(
    tdb$set_query(dims = list(Dept = c("A", "B")))
  )

  # verify selected ranges were set
  expect_equal(
    tiledb::selected_ranges(tdb$object),
    list(Dept = cbind(c("A", "B"), c("A", "B")))
  )

  # query result includes only selected ranges
  expect_equal(
    unique(tdb$object[]$Dept),
    c("A", "B")
  )

  # set attribute filter
  tdb <- TileDBArray$new(uri = uri, internal_use_only = "allowed_use")
  tdb$set_query(attr_filter = Admit == "Admitted")
  expect_true(all(tdb$object[]$Admit == "Admitted"))

  # update attribute filter
  tdb$set_query(attr_filter = Admit != "Admitted")
  expect_true(all(tdb$object[]$Admit == "Rejected"))

  # reset attribute filter
  tdb$reset_query()
  expect_length(tdb$object[]$Admit, nrow(df))

  ## shape
  expect_equal(tdb$ndim(), 2)

})
