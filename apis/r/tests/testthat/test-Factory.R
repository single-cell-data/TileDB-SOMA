test_that("DataFrame Factory", {
    uri <- tempfile()

    # check that straight use of new() errors, but 'with handshake' passes
    expect_error(SOMADataFrame$new(uri))
    expect_silent(d1 <- SOMADataFrame$new(uri, internal_use_only = "allowed_use"))

    # check creation of a DF
    asch <- create_arrow_schema()
    expect_error(d2 <- SOMADataFrameCreate(uri, asch)) # misses ind col name
    expect_error(d2 <- SOMADataFrameCreate(uri, index_column_names = "foo")) # misses scheme
    expect_silent(d2 <- SOMADataFrameCreate(uri, schema = asch, index_column_names = "foo"))
    tbl <- arrow::arrow_table(foo = 1L:10L, soma_joinid = 1L:10L, bar = sqrt(1:10),
                              baz = letters[1:10], schema = asch)
    d2$write(tbl)

    # check opening to read
    expect_silent(d3 <- SOMADataFrameOpen(uri))
    expect_silent(chk <- d3$read())
    expect_equal(tbl, chk)
})
