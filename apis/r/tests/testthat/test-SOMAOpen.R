test_that("SOMAOpen", {
    skip_if_not_installed("pbmc3k.tiledb")      # a Suggests: pre-package 3k PBMC data

    tdir <- tempfile()
    tgzfile <- system.file("raw-data", "soco-pbmc3k.tar.gz", package="pbmc3k.tiledb")
    untar(tarfile = tgzfile, exdir = tdir)
    uri <- file.path(tdir, "soco")

    expect_error(SOMAOpen(tdir)) # we cannot open directories are neither TileDB Array nor Group

    expect_true(inherits(SOMAOpen(uri), "SOMACollection"))

    expect_true(inherits(SOMAOpen(file.path(uri, "pbmc3k_processed")), "SOMAExperiment"))

    expect_true(inherits(SOMAOpen(file.path(uri, "pbmc3k_processed", "ms")), "SOMACollection"))
    expect_true(inherits(SOMAOpen(file.path(uri, "pbmc3k_processed", "obs")), "SOMADataFrame"))

    expect_true(inherits(SOMAOpen(file.path(uri, "pbmc3k_processed", "ms", "raw")), "SOMAMeasurement"))
    expect_true(inherits(SOMAOpen(file.path(uri, "pbmc3k_processed", "ms", "RNA")), "SOMAMeasurement"))

    expect_true(inherits(SOMAOpen(file.path(uri, "pbmc3k_processed", "ms", "RNA", "obsm", "X_draw_graph_fr")), "SOMASparseNDArray"))

    expect_true(inherits(SOMAOpen(file.path(uri, "pbmc3k_processed", "ms", "raw", "X", "data")), "SOMASparseNDArray"))

})
