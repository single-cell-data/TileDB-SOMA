test_that("Stats generation", {
    uri <- tempfile()
    create_and_populate_soma_dataframe(uri)

    sdf <- SOMADataFrame$new(uri)
    tiledb_stats_enable()
    arr <- sdf$read()
    txt <- tiledb_stats_dump()
    expect_true(nchar(txt) > 1000) # cannot parse JSON without a JSON package

    tiledb_stats_reset()
    txt <- tiledb_stats_dump()
    expect_true(nchar(txt) < 100) # almost empty JSON string

})
