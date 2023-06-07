test_that("Arrow Interface from SOMAArrayReader", {
    library(arrow)
    library(tiledb)

    uri <- extract_dataset("soma-dataframe-pbmc3k-processed-obs")

    columns <- c("n_counts", "n_genes", "louvain")

    z <- soma_array_reader(uri, columns)
    tb <- soma_array_to_arrow_table(z)
    expect_true(is_arrow_table(tb))
    rb <- arrow::as_record_batch(tb)
    expect_true(is_arrow_record_batch(rb))

    tb1 <- soma_array_to_arrow_table(soma_array_reader(uri, columns))
    expect_equal(tb1$num_rows, 2638)

    arr <- tiledb_array(uri)         # need array for schema access to qc parser
    qc <- parse_query_condition(n_counts < 1000 && n_genes >= 400, ta = arr)
    tb2 <- soma_array_to_arrow_table(soma_array_reader(uri, columns, qc@ptr))

    expect_equal(tb2$num_rows, 47)
    expect_true(all(tb2$n_counts < 1000))
    expect_true(all(tb2$n_genes >= 400))


    # read everything
    tb3 <- soma_array_to_arrow_table(soma_array_reader(uri))

    expect_equal(tb3$num_rows, 2638)
    expect_equal(tb3$num_columns, 6)

    # read a subset of rows and columns
    tb4 <- soma_array_reader(uri = uri,
                colnames = c("obs_id", "percent_mito", "n_counts", "louvain"),
                dim_ranges = list(soma_joinid = rbind(bit64::as.integer64(c(1000, 1004)),
                                                  bit64::as.integer64(c(2000, 2004)))),
                dim_points=list(soma_joinid = bit64::as.integer64(seq(0, 100, by = 20)))) |>
        soma_array_to_arrow_table()

    expect_equal(tb4$num_rows, 16)
    expect_equal(tb4$num_columns, 4)
})


test_that("SOMAArrayReader result order", {
    uri <- withr::local_tempdir("soma-dense-ndarray")
    ndarray <- SOMADenseNDArrayCreate(uri, arrow::int32(), shape = c(4, 4))

    M <- matrix(1:16, 4, 4)
    ndarray$write(M)
    ndarray$close()

    M1 <- soma_array_reader(uri = uri, result_order = "auto") |>
        soma_array_to_arrow_table()
    expect_equal(M, matrix(M1$soma_data, 4, 4, byrow = TRUE))

    M2 <- soma_array_reader(uri = uri, result_order = "row-major") |>
        soma_array_to_arrow_table()
    expect_equal(M, matrix(M2$soma_data, 4, 4, byrow = TRUE))

    M3 <- soma_array_reader(uri = uri, result_order = "column-major") |>
        soma_array_to_arrow_table()
    expect_equal(M, matrix(M3$soma_data, 4, 4, byrow = FALSE))
})
