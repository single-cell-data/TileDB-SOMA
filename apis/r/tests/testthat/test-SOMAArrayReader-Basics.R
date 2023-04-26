test_that("Basic SOMAArrayReader", {
    uri <- extract_dataset("soma-dataframe-pbmc3k-processed-obs")

    df <- arrow_to_dt(soma_array_reader(uri))
    expect_equal(nrow(df), 2638L)
    expect_equal(ncol(df), 6L)

    columns <- c("n_counts", "n_genes", "louvain")
    z <- soma_array_reader(uri, columns)
    expect_true(inherits(z, "list"))
})

test_that("Errors on large data, passes with budget", {
    skip_if_not_installed("pbmc3k.tiledb")      # a Suggests: pre-package 3k PBMC data
                                                # see https://ghrr.github.io/drat/

    tdir <- tempfile()
    tgzfile <- system.file("raw-data", "soco-pbmc3k.tar.gz", package="pbmc3k.tiledb")
    untar(tarfile = tgzfile, exdir = tdir)

    uri <- file.path(tdir, "soco", "pbmc3k_processed", "ms", "RNA", "X", "data")
    expect_true(dir.exists(uri))

    ## error under 'normal' read with 16mb default budget
    sdf1 <- SOMADataFrameOpen(uri)
    expect_error(sdf1$read())

    ## pass with budget of 45mb
    ctx <- tiledbsoma::SOMATileDBContext$new(c(soma.init_buffer_bytes="45000000"))
    sdf2 <- SOMADataFrameOpen(uri, tiledbsoma_ctx = ctx)
    expect_silent(sdf2$read())
})
