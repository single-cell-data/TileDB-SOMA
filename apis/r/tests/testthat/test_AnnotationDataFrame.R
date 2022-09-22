test_that("annotation dataframe can be stored and retrieved", {
  uri <- withr::local_tempdir("annot-df")

  df <- data.frame(
    character = c("A", "B", "C"),
    double = c(1.0, 2.0, 3.0),
    integer = c(1L, 2L, 3L),
    logical = c(TRUE, FALSE, TRUE)
  )

  annotdf <- AnnotationDataframe$new(uri)
  expect_true(inherits(annotdf, "AnnotationDataframe"))
  expect_error(
    annotdf$from_dataframe(df, index_col = "index"),
    "'x' must have character row names"
  )

  rownames(df) <- c("a", "b", "c")
  annotdf$from_dataframe(df, index_col = "index")

  expect_true(dir.exists(annotdf$uri))
  expect_true(annotdf$exists())
  expect_s4_class(annotdf$tiledb_array(), "tiledb_array")
  expect_is(annotdf$object, "tiledb_array")

  # data types
  expect_equal(tiledb::datatype(annotdf$attributes()[["character"]]), "ASCII")
  expect_equal(tiledb::datatype(annotdf$attributes()[["double"]]), "FLOAT64")
  expect_equal(tiledb::datatype(annotdf$attributes()[["integer"]]), "INT32")
  expect_equal(tiledb::datatype(annotdf$attributes()[["logical"]]), "BOOL")

  # helpers
  expect_setequal(annotdf$ids(), rownames(df))

  # retrieved dataframe
  df2 <- annotdf$to_dataframe()
  expect_setequal(rownames(df2), rownames(df))
  expect_setequal(colnames(df2), colnames(df))
  expect_identical(df2, df)
})

test_that("an empty dataframe can be stored and retrieved", {
  uri <- withr::local_tempdir("annot-df-empty")
  df <- data.frame(row.names = letters)
  expect_length(df, 0)

  annotdf <- AnnotationDataframe$new(uri)
  annotdf$from_dataframe(df, index_col = "index")

  df2 <- annotdf$to_dataframe()
  expect_length(df2, 0)
  expect_setequal(rownames(df2), rownames(df))
})

test_that("updates overwrite existing cells", {
  uri <- withr::local_tempdir("annot-df-updates")
  df <- data.frame(row.names = "a", value = 1)

  annotdf <- AnnotationDataframe$new(uri)
  annotdf$from_dataframe(df, index_col = "index")
  expect_identical(annotdf$to_dataframe(), df)

  df2 <- data.frame(row.names = "a", value = 2)
  annotdf$from_dataframe(df2, index_col = "index")
  expect_identical(annotdf$to_dataframe(), df2)
})
