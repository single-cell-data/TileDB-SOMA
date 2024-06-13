test_that("DataFrame Factory", {
    skip_if(!extended_tests())
    uri <- tempfile()

    # Check that straight use of new() errors, but 'with handshake' passes
    expect_error(SOMADataFrame$new(uri))
    expect_silent(d1 <- SOMADataFrame$new(uri, internal_use_only = "allowed_use"))

    # Check creation of a DF
    asch <- create_arrow_schema(foo_first=FALSE)
    expect_silent(d2 <- SOMADataFrameCreate(uri, schema = asch))
    tbl <- arrow::arrow_table(soma_joinid = 1L:10L, foo = 1L:10L, bar = sqrt(1:10),
                              baz = letters[1:10], schema = asch)
    d2$write(tbl)

    # Check opening to read
    expect_silent(d3 <- SOMADataFrameOpen(uri))
    expect_silent(chk <- d3$read()$concat())
    expect_equal(tibble::as_tibble(tbl), tibble::as_tibble(chk))
})

test_that("DataFrame Factory with specified index_column_names", {
    skip_if(!extended_tests())
    uri <- tempfile()

    # Check creation of a DF
    asch <- create_arrow_schema()
    expect_error(d2 <- SOMADataFrameCreate(uri, index_column_names = "foo")) # misses schema
    expect_silent(d2 <- SOMADataFrameCreate(uri, schema = asch, index_column_names = "foo"))
    tbl <- arrow::arrow_table(foo = 1L:10L, soma_joinid = 1L:10L, bar = sqrt(1:10),
                              baz = letters[1:10], schema = asch)
    d2$write(tbl)

    # Check opening to read
    expect_silent(d3 <- SOMADataFrameOpen(uri))
    expect_equal(d3$mode(), "READ")
    expect_silent(chk <- d3$read()$concat())
    expect_equal(tibble::as_tibble(tbl), tibble::as_tibble(chk))
    d3$close()
    expect_equal(d3$mode(), "CLOSED")
})

test_that("SparseNDArray Factory", {
    skip_if(TRUE)
    skip_if(!extended_tests())
    uri <- tempfile()

    # check that straight use of new() errors, but 'with handshake' passes
    expect_error(SOMASparseNDArray$new(uri))
    expect_silent(s1 <- SOMASparseNDArray$new(uri, internal_use_only = "allowed_use"))

    # check creation of a sparse array
    expect_error(s2 <- SOMASparseNDArrayCreate(uri, arrow::int32())) # misses shape
    expect_error(s2 <- SOMASparseNDArrayCreate(uri, shape = c(10,10))) # misses type
    expect_silent(s2 <- SOMASparseNDArrayCreate(uri, arrow::int32(), shape = c(10,10)))
    mat <- create_sparse_matrix_with_int_dims(10, 10)
    s2$write(mat)

    # check opening to read
    expect_silent(s3 <- SOMASparseNDArrayOpen(uri))
    expect_equal(s3$mode(), "READ")

    #TODO test when sr_setup has an argument "result_order"
    #expect_silent(chk <- s3$read(result_order = "COL_MAJOR")$tables()$concat())
    #expect_identical(
    #    as.numeric(chk$GetColumnByName("soma_data")),
    #    ## need to convert to Csparsematrix first to get x values sorted appropriately
    #    as.numeric(as(mat, "CsparseMatrix")@x)
    #)
    s3$close()
    expect_equal(s3$mode(), "CLOSED")
})

test_that("DenseNDArray Factory", {
    skip_if(!extended_tests())
    uri <- tempfile()

    # check that straight use of new() errors, but 'with handshake' passes
    expect_error(SOMADenseNDArray$new(uri))
    expect_silent(s1 <- SOMADenseNDArray$new(uri, internal_use_only = "allowed_use"))

    # check creation of a sparse array
    expect_error(s2 <- SOMADenseNDArrayCreate(uri, arrow::int32())) # misses shape
    expect_error(s2 <- SOMADenseNDArrayCreate(uri, shape = c(10,10))) # misses type
    expect_silent(s2 <- SOMADenseNDArrayCreate(uri, arrow::int32(), shape = c(10,10)))
    mat <- create_dense_matrix_with_int_dims(10, 10)
    s2$write(mat)

    # check opening to read
    expect_silent(s3 <- SOMADenseNDArrayOpen(uri))
    expect_equal(s3$mode(), "READ")
    expect_silent(chk <- s3$read_dense_matrix())
    expect_equal(mat, chk)
    s3$close()
    expect_equal(s3$mode(), "CLOSED")
})

test_that("Collection Factory", {
    skip_if(!extended_tests())
    uri <- tempfile()

    # check that straight use of new() errors, but 'with handshake' passes
    expect_error(SOMACollection$new(uri))
    expect_silent(s1 <- SOMACollection$new(uri, internal_use_only = "allowed_use"))

    # check creation of a sparse array
    expect_silent(s2 <- SOMACollectionCreate(uri))

    # check opening to read
    expect_silent(s3 <- SOMACollectionOpen(uri))
    expect_equal(s3$mode(), "READ")
    s3$close()
    expect_equal(s3$mode(), "CLOSED")
})

test_that("Measurement Factory", {
    skip_if(!extended_tests())
    uri <- tempfile()

    # check that straight use of new() errors, but 'with handshake' passes
    expect_error(SOMAMeasurement$new(uri))
    expect_silent(s1 <- SOMAMeasurement$new(uri, internal_use_only = "allowed_use"))

    # check creation of a sparse array
    expect_silent(s2 <- SOMAMeasurementCreate(uri))

    # check opening to read
    expect_silent(s3 <- SOMAMeasurementOpen(uri))
    expect_equal(s3$mode(), "READ")
    s3$close()
    expect_equal(s3$mode(), "CLOSED")
})

test_that("Experiment Factory", {
    skip_if(!extended_tests())
    uri <- tempfile()

    # check that straight use of new() errors, but 'with handshake' passes
    expect_error(SOMAExperiment$new(uri))
    expect_silent(s1 <- SOMAExperiment$new(uri, internal_use_only = "allowed_use"))

    # check creation of a sparse array
    expect_silent(s2 <- SOMAExperimentCreate(uri))

    # check opening to read
    expect_silent(s3 <- SOMAExperimentOpen(uri))
    expect_equal(s3$mode(), "READ")
    s3$close()
    expect_equal(s3$mode(), "CLOSED")
})
