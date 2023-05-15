test_that("Stats generation", {
    uri <- tempfile()
    sdf <- create_and_populate_soma_dataframe(uri, mode = "READ")
    on.exit(sdf$close())

    tiledbsoma_stats_enable()
    arr <- sdf$read()
    txt <- tiledbsoma_stats_dump()
    expect_true(nchar(txt) > 1000) # cannot parse JSON without a JSON package

    tiledbsoma_stats_reset()
    txt <- tiledbsoma_stats_dump()
    expect_true(nchar(txt) < 100) # almost empty JSON string

})
