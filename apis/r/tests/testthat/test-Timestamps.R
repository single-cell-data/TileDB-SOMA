test_that("SOMADataFrame", {
    uri <- tempfile()

    sch <- arrow::schema(arrow::field("soma_joinid", arrow::int64()),
                         arrow::field("int", arrow::int32()),
                         arrow::field("str", arrow::dictionary(index_type = arrow::int8(),
                                                               value_type = arrow::utf8())))

    ## create at t = 1
    ts1 <- as.POSIXct(1, tz = "UTC", origin = "1970-01-01")
    sdf <- tiledbsoma::SOMADataFrameCreate(uri, sch, tiledb_timestamp = ts1)

    ## write part1 at t = 2
    dat2 <- arrow::arrow_table(soma_joinid = bit64::as.integer64(1L:5L),
                               int = 101:105L,
                               str = factor(c('a','b','b','a','b')))
    ts2 <- as.POSIXct(2, tz = "UTC", origin = "1970-01-01")
    expect_silent(sdf$write(dat2, ts2))

    ## write part2 at t = 3
    dat3 <- arrow::arrow_table(soma_joinid = bit64::as.integer64(6L:10L),
                               int = 106:110L,
                               str = factor(c('c','b','c','c','b')))
    ts3 <- as.POSIXct(3, tz = "UTC", origin = "1970-01-01")
    expect_silent(sdf$write(dat3, ts3))
    sdf$close()

    ## read all
    sdf <- tiledbsoma::SOMADataFrameOpen(uri)
    res <- tibble::as_tibble(sdf$read()$concat())
    expect_equal(dim(res), c(10, 3)) 					# two writes lead to 10x3 data
    expect_equal(levels(res$str), c("a", "b", "c"))     # string variable has three values

    ## read before data is written (tiledb_timestamp_range = <origin, ts1>)
    sdf <- tiledbsoma::SOMADataFrameOpen(uri, tiledb_timestamp = ts1)
    res <- sdf$read()$concat()
    expect_equal(dim(res), c(0, 3))

    ## read at ts2 (tiledb_timestamp_range = <origin, ts2>)
    sdf <- tiledbsoma::SOMADataFrameOpen(uri, tiledb_timestamp = ts2)
    res <- tibble::as_tibble(sdf$read()$concat())
    expect_equal(dim(res), c(5, 3))
    expect_equal(max(res$int), 105L)
    expect_equal(range(res$int), c(101L,105L))

    ## read at ts3 (tiledb_timestamp_range = <origin, ts3>)
    sdf <- tiledbsoma::SOMADataFrameOpen(uri, tiledb_timestamp = ts3)
    res <- tibble::as_tibble(sdf$read()$concat())
    expect_equal(dim(res), c(10L, 3L))
    expect_equal(max(res$int), 110L)
    expect_equal(range(res$int), c(101L, 110L))

    ## read after ts3 (tiledb_timestamp_range = <origin, ts3 + 1>)
    sdf <- tiledbsoma::SOMADataFrameOpen(uri, tiledb_timestamp = ts3 + 1)
    res <- tibble::as_tibble(sdf$read()$concat())
    res <- sdf$read()$concat()
    expect_equal(dim(res), c(10L, 3L))
})

test_that("SOMANDSparseArray", {
    uri <- tempfile()

    ## create at t = 2, also writes at t = 2 as timestamp is cached
    ts1 <- as.POSIXct(2, tz = "UTC", origin = "1970-01-01")
    snda <- SOMASparseNDArrayCreate(
        uri,
        arrow::int32(),
        shape = c(10, 10),
        tiledb_timestamp = ts1
    )
    mat <- create_sparse_matrix_with_int_dims(10, 10)
    snda$write(mat) 		# write happens at create time

    ## read at t = 3, expect all rows as read is from (0, 3)
    ts2 <- as.POSIXct(3, tz = "UTC", origin = "1970-01-01")
    snda <- tiledbsoma::SOMASparseNDArrayOpen(uri, tiledb_timestamp = ts2)
    res <- snda$read()$tables()$concat()
    expect_equal(dim(res), c(60,3))

    ## read at t = 1, expect zero rows as read is from (0, 1)
    ts3 <- as.POSIXct(1, tz = "UTC", origin = "1970-01-01")
    snda <- tiledbsoma::SOMASparseNDArrayOpen(uri, tiledb_timestamp = ts3)
    res <- snda$read()$tables()$concat()
    expect_equal(dim(res), c(0,3))

})

test_that("SOMANDDenseArray", {
    uri <- tempfile()

    ## create at t = 2, also writes at t = 2 as timestamp is cached
    ts1 <- as.POSIXct(2, tz = "UTC", origin = "1970-01-01")
    dnda <- SOMADenseNDArrayCreate(
        uri,
        type = arrow::int32(),
        shape = c(5, 2),
        tiledb_timestamp = ts1
    )
    mat <- create_dense_matrix_with_int_dims(5, 2)
    expect_silent(dnda$write(mat)) 		# write happens at create time

    ## read at t = 3, expect all rows as read is from (0, 3)
    ts2 <- as.POSIXct(3, tz = "UTC", origin = "1970-01-01")
    dnda <- tiledbsoma::SOMADenseNDArrayOpen(uri, tiledb_timestamp = ts2)
    res <- dnda$read_arrow_table()
    expect_equal(dim(res), c(10,1))

    ## read at t = 1, expect zero rows as read is from (0, 1)
    ## NB that this requires a) nullable arrow schema and b) na.omit() on result
    ## It also works without a) as the fill value is also NA
    ts3 <- as.POSIXct(1, tz = "UTC", origin = "1970-01-01")
    dnda <- tiledbsoma::SOMADenseNDArrayOpen(uri, tiledb_timestamp = ts3)
    res <- dnda$read_arrow_table()
    #print(res$soma_data)
    vec <- tibble::as_tibble(res)
    expect_equal(dim(na.omit(vec)), c(0,1))

})
