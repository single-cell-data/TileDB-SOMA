test_that("SOMADataFrame", {
    uri <- tempfile()

    sch <- arrow::schema(arrow::field("soma_joinid", arrow::int64()),
                         arrow::field("int", arrow::int32()),
                         arrow::field("str", arrow::dictionary(index_type = arrow::int8(),
                                                               value_type = arrow::utf8())))

    ## create at t = 1
    ts1 <- rep(as.POSIXct(1, tz="UTC"), 2)
    sdf <- tiledbsoma::SOMADataFrameCreate(uri, sch, tiledb_timestamp=ts1[1])

    ## write part1 at t = 2
    dat2 <- arrow::arrow_table(soma_joinid = bit64::as.integer64(1L:5L),
                               int = 101:105L,
                               str = factor(c('a','b','b','a','b')))
    ts2 <- rep(as.POSIXct(2, tz="UTC"), 2)
    expect_silent(sdf$write(dat2, ts2))

    ## write part2 at t = 3
    dat3 <- arrow::arrow_table(soma_joinid = bit64::as.integer64(6L:10L),
                               int = 106:110L,
                               str = factor(c('c','b','c','c','b')))
    ts3 <- rep(as.POSIXct(3, tz="UTC"), 2)
    expect_silent(sdf$write(dat3, ts3))
    sdf$close()

    res <- tiledb::tiledb_array(uri, return="data.frame")[]
    expect_equal(dim(res), c(10, 3)) 					# two writes lead to 10x3 data
    expect_equal(levels(res$str), c("a", "b", "c"))     # factor variable has three levels
})
