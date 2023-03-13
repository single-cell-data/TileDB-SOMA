test_that("read and write array without legacy mode", {
  on.exit(tiledb::tiledb_ctx(tiledb::tiledb_config()))
  uri <- withr::local_tempdir("df-no-legacy")

  expect_null(tiledb_ctx_get_key(TILEDB_LEGACY_KEY))

  df <- data.frame(
    character = c("A", "B", "C", NA_character_),
    row.names = c("a", "b", "c", "d")
  )

  annotdf <- AnnotationDataframe$new(uri)
  annotdf$from_dataframe(df, index_col = "index")

  # verify legacy metadata tag is present and reads are valid
  expect_equal(annotdf$get_metadata(SOMA_LEGACY_VALIDITY_KEY), "false")
  expect_equal(annotdf$to_dataframe()$character, df$character)

  # config param is still unset
  expect_null(tiledb_ctx_get_key(TILEDB_LEGACY_KEY))

  # enabling legacy throws error on reads
  tiledb_ctx_set_key(TILEDB_LEGACY_KEY, "true")
  expect_error(
    annotdf$to_dataframe(),
    "Legacy mode is enabled but this array was created without it"
  )

  # enabling legacy mode throws error on init
  expect_error(
    AnnotationDataframe$new(uri),
    "Legacy mode is enabled but this array was created without it"
  )
})
