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

test_that("SparseNDArray Factory", {
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
    expect_silent(chk <- s3$read_arrow_table(result_order = "COL_MAJOR"))
    expect_identical(
        as.numeric(chk$GetColumnByName("soma_data")),
        ## need to convert to Csparsematrix first to get x values sorted appropriately
        as.numeric(as(mat, "CsparseMatrix")@x)
    )

})

test_that("SparseNDArray Factory", {
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
    expect_silent(chk <- s3$read_dense_matrix())
    expect_equal(mat, chk)
})

test_that("Collection Factory", {
    uri <- tempfile()

    # check that straight use of new() errors, but 'with handshake' passes
    expect_error(SOMACollection$new(uri))
    expect_silent(s1 <- SOMACollection$new(uri, internal_use_only = "allowed_use"))

    # check creation of a sparse array
    expect_silent(s2 <- SOMACollectionCreate(uri))

    # check opening to read
    expect_silent(s3 <- SOMACollectionOpen(uri))
})

test_that("Measurement Factory", {
    uri <- tempfile()

    # check that straight use of new() errors, but 'with handshake' passes
    expect_error(SOMAMeasurement$new(uri))
    expect_silent(s1 <- SOMAMeasurement$new(uri, internal_use_only = "allowed_use"))

    # check creation of a sparse array
    expect_silent(s2 <- SOMAMeasurementCreate(uri))

    # check opening to read
    expect_silent(s3 <- SOMAMeasurementOpen(uri))
})

test_that("Experiment Factory", {
    uri <- tempfile()

    # check that straight use of new() errors, but 'with handshake' passes
    expect_error(SOMAExperiment$new(uri))
    expect_silent(s1 <- SOMAExperiment$new(uri, internal_use_only = "allowed_use"))

    # check creation of a sparse array
    expect_silent(s2 <- SOMAExperimentCreate(uri))

    # check opening to read
    expect_silent(s3 <- SOMAExperimentOpen(uri))
})
