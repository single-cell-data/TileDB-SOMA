test_that("Blockwise iterator for arrow tables", {
    skip_if(!extended_tests() || covr_tests())
    skip_if_not_installed("pbmc3k.tiledb")      # a Suggests: pre-package 3k PBMC data
                                                # see https://ghrr.github.io/drat/

    tdir <- tempfile()
    tgzfile <- system.file("raw-data", "soco-pbmc3k.tar.gz", package="pbmc3k.tiledb")
    untar(tarfile = tgzfile, exdir = tdir)

    uri <- file.path(tdir, "soco", "pbmc3k_processed")
    expect_true(dir.exists(uri))

    ax <- 0
    sz <- 1000L
    expqry <- SOMAExperimentOpen(uri)
    axqry <- expqry$axis_query("RNA")
    xrqry <- axqry$X("data")

    expect_error(xrqry$blockwise(axis=2))
    expect_error(xrqry$blockwise(size=-100))

    expect_s3_class(
        bi <- xrqry$blockwise(axis=ax, size=sz),
        "SOMASparseNDArrayBlockwiseRead"
    )

    expect_s3_class(it <- bi$tables(), "BlockwiseTableReadIter")
    expect_false(it$read_complete())

    for (i in seq.int(1L, ceiling(it$coords_axis$length() / it$coords_axis$stride))) {
        at <- it$read_next()
        expect_s3_class(at, "ArrowTabular")
    }
    expect_true(it$read_complete())

    rm(bi, it, xrqry, axqry)
    axqry <- expqry$axis_query("RNA")
    xrqry <- axqry$X("data")
    bi <- xrqry$blockwise(axis=ax, size=sz)
    it <- bi$tables()
    at <- it$concat()
    expect_s3_class(at, "Table")
    expect_s3_class(at, "ArrowTabular")
    expect_equal(dim(at), c(4848644, 3))
})

test_that("Blockwise iterator for sparse matrices", {
    skip_if(!extended_tests() || covr_tests())
    skip_if_not_installed("pbmc3k.tiledb")      # a Suggests: pre-package 3k PBMC data
                                                # see https://ghrr.github.io/drat/

    tdir <- tempfile()
    tgzfile <- system.file("raw-data", "soco-pbmc3k.tar.gz", package="pbmc3k.tiledb")
    untar(tarfile = tgzfile, exdir = tdir)

    uri <- file.path(tdir, "soco", "pbmc3k_processed")
    expect_true(dir.exists(uri))

    ax <- 0
    sz <- 1000L
    expqry <- SOMAExperimentOpen(uri)
    axqry <- expqry$axis_query("RNA")
    xrqry <- axqry$X("data")

    expect_error(xrqry$blockwise(axis=2))
    expect_error(xrqry$blockwise(size=-100))

    expect_s3_class(
        bi <- xrqry$blockwise(axis=ax, size=sz),
        "SOMASparseNDArrayBlockwiseRead"
    )

    expect_s3_class(it <- bi$sparse_matrix(), "BlockwiseSparseReadIter")
    expect_false(it$read_complete())

    for (i in seq.int(1L, ceiling(it$coords_axis$length() / it$coords_axis$stride))) {
        at <- it$read_next()
        expect_s4_class(at, "dgTMatrix")
    }
    expect_true(it$read_complete())

    rm(bi, it, xrqry, axqry)
    axqry <- expqry$axis_query("RNA")
    xrqry <- axqry$X("data")
    bi <- xrqry$blockwise(axis=ax, size=sz)
    it <- bi$sparse_matrix()
    at <- it$concat()
    expect_s4_class(at, "dgTMatrix")
    expect_equal(dim(at), c(2638, 1838))
})
