test_that("Stats generation", {
    uri <- tempfile()
    create_and_populate_soma_dataframe(uri)

    sdf <- SOMADataFrame$new(uri)
    sdf$stats_enable()
    arr <- sdf$read()
    txt <- sdf$stats_dump()
    expect_true(nchar(txt) > 1000) # cannot parse JSON without a JSON package

    sdf$stats_reset()
    txt <- sdf$stats_dump()
    expect_true(nchar(txt) < 100) # almost empty JSON string

})
