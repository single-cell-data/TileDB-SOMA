test_that("Arrow Interface from SOMAArrayReader", {
    library(arrow)
    library(tiledb)

    skip_if_not_installed("dplyr")      # a Suggests
    uri <- extract_dataset("soma-dataframe-pbmc3k-processed-obs")

    columns <- c("n_counts", "n_genes", "louvain")

    z <- soma_array_reader(uri, columns)
    tb <- to_arrow_table(z)
    expect_true(inherits(tb, "Table"))
    rb <- arrow::as_record_batch(tb)  #arch::from_arch_array(z, arrow::RecordBatch)
    expect_true(inherits(rb, "RecordBatch"))


    soma_array_reader(uri, columns) |>
        to_arrow_table() |>
        dplyr::collect() -> D
    expect_equal(nrow(D), 2638)

    arr <- tiledb_array(uri)                # need array for schema access to qc parser
    qc <- parse_query_condition(n_counts < 1000 && n_genes >= 400, ta=arr)
    soma_array_reader(uri, columns, qc@ptr) |>
        to_arrow_table() |>
        dplyr::collect() -> D

    expect_equal(nrow(D), 47)
    expect_true(all(D$n_counts < 1000))
    expect_true(all(D$n_genes >= 400))


    soma_array_reader(uri) |>              # read everything
        to_arrow_table() |>
        dplyr::collect() -> D
    expect_equal(nrow(D), 2638)
    expect_equal(ncol(D), 6)

    soma_array_reader(uri = uri,
                colnames = c("obs_id", "percent_mito", "n_counts", "louvain"),
                dim_ranges=list(soma_joinid=rbind(bit64::as.integer64(c(1000, 1004)),
                                                  bit64::as.integer64(c(2000, 2004)))),
                dim_points=list(soma_joinid=bit64::as.integer64(seq(0, 100, by=20)))) |>
        to_arrow_table() |>
        dplyr::collect() -> D
    expect_equal(nrow(D), 16)
    expect_equal(ncol(D), 4)


    uri <- tempfile()
    ndarray <- SOMADenseNDArrayCreate(uri, arrow::int32(), shape = c(4, 4))
    M <- matrix(1:16, 4, 4)
    ndarray$write(M)
    ndarray$close()

    M1 <- soma_array_reader(uri = uri, result_order = "auto") |>
        to_arrow_table() |>
        dplyr::collect()
    expect_equal(M, matrix(M1$soma_data, 4, 4, byrow=TRUE))

    M2 <- soma_array_reader(uri = uri, result_order = "row-major") |>
        to_arrow_table() |>
        dplyr::collect()
    expect_equal(M, matrix(M2$soma_data, 4, 4, byrow=TRUE))

    M3 <- soma_array_reader(uri = uri, result_order = "column-major") |>
        to_arrow_table() |>
        dplyr::collect()
    expect_equal(M, matrix(M3$soma_data, 4, 4, byrow=FALSE))


})
