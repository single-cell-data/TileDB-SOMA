test_that("SOMADataFrame", {
    uri <- tempfile()

    sch <- arrow::schema(arrow::field("soma_joinid", arrow::int64()),
                         arrow::field("int", arrow::int32()),
                         #arrow::field("str", arrow::dictionary(index_type = arrow::int8(),
                         #                                      value_type = arrow::utf8())))
                         arrow::field("str", arrow::utf8()))

    ## create at t = 1
    ts1 <- rep(as.POSIXct(1, tz="UTC"), 2)
    sdf <- tiledbsoma::SOMADataFrameCreate(uri, sch, tiledb_timestamp=ts1[1])

    ## write part1 at t = 2
    dat2 <- arrow::arrow_table(soma_joinid = bit64::as.integer64(1L:5L),
                               int = 101:105L,
                               #str = factor(c('a','b','b','a','b')))
                               str = c('a','b','b','a','b'))
    ts2 <- rep(as.POSIXct(2, tz="UTC"), 2)
    expect_silent(sdf$write(dat2, ts2))

    ## write part2 at t = 3
    dat3 <- arrow::arrow_table(soma_joinid = bit64::as.integer64(6L:10L),
                               int = 106:110L,
                               #str = factor(c('c','b','c','c','b')))
                               str = c('c','b','c','c','b'))
    ts3 <- rep(as.POSIXct(3, tz="UTC"), 2)
    expect_silent(sdf$write(dat3, ts3))
    sdf$close()

    ## read all
    sdf <- tiledbsoma::SOMADataFrameOpen(uri)
    res <- tibble::as_tibble(sdf$read()$concat())
    expect_equal(dim(res), c(10, 3)) 					# two writes lead to 10x3 data
    expect_equal(unique(res$str), c("a", "b", "c"))     # string variable has three values

    ## read before data is written
    sdf <- tiledbsoma::SOMADataFrameOpen(uri, tiledb_timestamp=ts1)
    res <- sdf$read()$concat()
    expect_equal(dim(res), c(0, 3))

    ## read at ts2
    sdf <- tiledbsoma::SOMADataFrameOpen(uri, tiledb_timestamp=ts2)
    res <- tibble::as_tibble(sdf$read()$concat())
    expect_equal(dim(res), c(5, 3))
    expect_equal(max(res$int), 105L)
    expect_equal(range(res$int), c(101L,105L))

    ## read at ts3
    sdf <- tiledbsoma::SOMADataFrameOpen(uri, tiledb_timestamp=ts3)
    res <- tibble::as_tibble(sdf$read()$concat())
    expect_equal(dim(res), c(5, 3))
    expect_equal(max(res$int), 110L)
    expect_equal(range(res$int), c(106L,110L))

    ## read after ts3
    sdf <- tiledbsoma::SOMADataFrameOpen(uri, tiledb_timestamp=ts3 + c(1,1))
    res <- tibble::as_tibble(sdf$read()$concat())
    res <- sdf$read()$concat()
    expect_equal(dim(res), c(0, 3))
})
