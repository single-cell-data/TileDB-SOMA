test_that("Arrow Interface from SOMAReader", {
    library(arch)
    library(arrow)
    library(tiledb)

    skip_if_not_installed("dplyr")      # a Suggests

    tdir <- tempfile()
    tgzfile <- system.file("raw-data", "soco-pbmc3k_processed-obs.tar.gz", package="tiledbsoma")
    untar(tarfile = tgzfile, exdir = tdir)

    uri <- file.path(tdir, "obs")
    columns <- c("n_counts", "n_genes", "louvain")

    z <- export_arrow_array(uri, columns)
    rb <- arch::from_arch_array(z, arrow::RecordBatch)
    expect_true(inherits(rb, "RecordBatch"))
    tb <- arrow::as_arrow_table(arch::from_arch_array(z, arrow::RecordBatch))
    expect_true(inherits(tb, "Table"))

    export_arrow_array(uri, columns, loglevel="warn") |>
        arch::from_arch_array(arrow::RecordBatch) |>
        arrow::as_arrow_table() |>
        dplyr::mutate(louvain = cast(louvain, arrow::string())) |>
        dplyr::collect() -> D
    expect_equal(nrow(D), 2638)

    arr <- tiledb_array(uri)                # need array for schema access to qc parser
    qc <- parse_query_condition(n_counts < 1000 && n_genes >= 400, ta=arr)
    export_arrow_array(uri, columns, qc@ptr, "warn") |>
        arch::from_arch_array(arrow::RecordBatch) |>
        arrow::as_arrow_table() |>
        mutate(louvain = cast(louvain, string())) |>
        collect() -> D

    expect_equal(nrow(D), 47)
    expect_true(all(D$n_counts < 1000))
    expect_true(all(D$n_genes >= 400))



export_arrow_array(uri) |>              # read everything
    arch::from_arch_array(arrow::RecordBatch) |>
    arrow::as_arrow_table() |>
    collect() -> D
expect_equal(nrow(D), 2638)
expect_equal(ncol(D), 6)

export_arrow_array(uri = uri,
                   colnames = c("obs_id", "percent_mito", "n_counts", "louvain"),
                   dim_ranges=list(soma_rowid=rbind(bit64::as.integer64(c(1000, 1004)),
                                                    bit64::as.integer64(c(2000, 2004)))),
                   dim_points=list(soma_rowid=bit64::as.integer64(seq(0, 100, by=20)))) |>
    arch::from_arch_array(arrow::RecordBatch) |>
    arrow::as_arrow_table() |>
    mutate(louvain = cast(louvain, string()),
           obs_id = cast(obs_id, string())) |>
    collect() -> D
expect_equal(nrow(D), 16)
expect_equal(ncol(D), 4)

})
