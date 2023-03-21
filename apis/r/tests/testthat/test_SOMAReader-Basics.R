test_that("Basic SOMAArrayReader", {
    tdir <- tempfile()
    tgzfile <- system.file("raw-data", "soco-pbmc3k_processed-obs.tar.gz", package="tiledbsoma")
    untar(tarfile = tgzfile, exdir = tdir)

    uri <- file.path(tdir, "obs")

    df <- arrow_to_dt(soma_reader(uri))
    expect_equal(nrow(df), 2638L)
    expect_equal(ncol(df), 6L)

    columns <- c("n_counts", "n_genes", "louvain")
    z <- soma_reader(uri, columns)
    expect_true(inherits(z, "list"))
})
