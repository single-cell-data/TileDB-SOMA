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

    bi <- xrqry$blockwise(axis=ax, size=sz)
    expect_true(inherits(bi, "SOMASparseNDArrayBlockwiseRead"))

    it <- bi$tables()
    expect_true(inherits(it, "BlockwiseTableReadIter"))
    expect_false(it$read_complete())

    at <- it$read_next()
    expect_true(inherits(at, "Table") && inherits(at, "ArrowTabular"))
    at <- it$read_next()
    expect_true(inherits(at, "Table") && inherits(at, "ArrowTabular"))
    at <- it$read_next()
    expect_true(inherits(at, "Table") && inherits(at, "ArrowTabular"))
    expect_true(it$read_complete())

    rm(bi, it, xrqry, axqry)
    axqry <- expqry$axis_query("RNA")
    xrqry <- axqry$X("data")
    bi <- xrqry$blockwise(axis=ax, size=sz)
    it <- bi$tables()
    at <- it$concat()
    expect_true(inherits(at, "Table") && inherits(at, "ArrowTabular"))
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

    bi <- xrqry$blockwise(axis=ax, size=sz)
    expect_true(inherits(bi, "SOMASparseNDArrayBlockwiseRead"))

    it <- bi$sparse_matrix()
    expect_true(inherits(it, "BlockwiseSparseReadIter"))
    expect_false(it$read_complete())

    at <- it$read_next()
    expect_true(inherits(at, "dgTMatrix"))
    at <- it$read_next()
    expect_true(inherits(at, "dgTMatrix"))
    at <- it$read_next()
    expect_true(inherits(at, "dgTMatrix"))
    expect_true(it$read_complete())

    rm(bi, it, xrqry, axqry)
    axqry <- expqry$axis_query("RNA")
    xrqry <- axqry$X("data")
    bi <- xrqry$blockwise(axis=ax, size=sz)
    it <- bi$sparse_matrix()
    at <- it$concat()
    expect_true(inherits(at, "dgTMatrix"))
    expect_equal(dim(at), c(2638, 1838))
})
