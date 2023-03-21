test_that("Iterated Interface from SOMAReader", {
    skip_if_not_installed("pbmc3k.tiledb")      # a Suggests: pre-package 3k PBMC data

    library(arrow)
    library(bit64)
    library(data.table)
    library(dplyr)
    library(tiledb)

    tdir <- tempfile()
    tgzfile <- system.file("raw-data", "soco-pbmc3k.tar.gz", package="pbmc3k.tiledb")
    untar(tarfile = tgzfile, exdir = tdir)

    uri <- file.path(tdir, "soco", "pbmc3k_processed", "ms", "RNA", "X", "data")
    expect_true(dir.exists(uri))

    ctx <- tiledb_ctx()
    sr <- sr_setup(uri, config=as.character(config(ctx)), loglevel="warn")
    expect_true(inherits(sr, "externalptr"))
    rl <- data.frame()
    while (!tiledbsoma:::sr_complete(sr)) {
        dat <- sr_next(sr)
        D <- tiledbsoma:::arrow_to_dt(dat)
        expect_true(nrow(D) > 0)
        expect_true(inherits(D, "data.table"))
        rl <- rbind(rl, D)
    }
    expect_true(inherits(rl, "data.table"))
    expect_equal(nrow(rl), 4848644)
    expect_equal(ncol(rl), 3)

    sr <- sr_setup(uri, config=as.character(config(ctx)), dim_points=list(soma_dim_0=as.integer64(1)))
    expect_true(inherits(sr, "externalptr"))
    rl <- data.frame()
    while (!tiledbsoma:::sr_complete(sr)) {
        dat <- sr_next(sr)
        D <- tiledbsoma:::arrow_to_dt(dat)
        expect_true(nrow(D) > 0)
        expect_true(inherits(D, "data.table"))
        rl <- rbind(rl, D)
    }
    expect_true(inherits(rl, "data.table"))
    expect_equal(nrow(rl), 1838)
    expect_equal(ncol(rl), 3)

    sr <- sr_setup(uri, config=as.character(config(ctx)), dim_range=list(soma_dim_1=cbind(as.integer64(1),as.integer64(2))))
    expect_true(inherits(sr, "externalptr"))
    rl <- data.frame()
    while (!tiledbsoma:::sr_complete(sr)) {
        dat <- sr_next(sr)
        D <- tiledbsoma:::arrow_to_dt(dat)
        expect_true(nrow(D) > 0)
        expect_true(inherits(D, "data.table"))
        rl <- rbind(rl, D)
    }
    expect_true(inherits(rl, "data.table"))
    expect_equal(nrow(rl), 5276)
    expect_equal(ncol(rl), 3)

    ## test completeness predicate on shorter data
    tdir <- tempfile()
    tgzfile <- system.file("raw-data", "soco-pbmc3k_processed-obs.tar.gz", package="tiledbsoma")
    untar(tarfile = tgzfile, exdir = tdir)
    uri <- file.path(tdir, "obs")
    sr <- sr_setup(uri, config=as.character(config(ctx)))

    expect_false(tiledbsoma:::sr_complete(sr))
    dat <- dat <- sr_next(sr)
    expect_true(tiledbsoma:::sr_complete(sr))
})

test_that("Iterated Interface from SOMA Classes", {
    skip_if_not_installed("pbmc3k.tiledb")      # a Suggests: pre-package 3k PBMC data

    tdir <- tempfile()
    tgzfile <- system.file("raw-data", "soco-pbmc3k.tar.gz", package="pbmc3k.tiledb")
    untar(tarfile = tgzfile, exdir = tdir)
    uri <- file.path(tdir, "soco", "pbmc3k_processed", "ms", "RNA", "X", "data")

    ## parameterize test
    test_cases <- c("data.frame", "sparse", "dense")

    for (tc in test_cases) {
        sdf <- switch(tc,
                      data.frame = SOMADataFrame$new(uri, internal_use_only = "allowed_use"),
                      sparse = SOMASparseNDArray$new(uri, internal_use_only = "allowed_use"),
                      dense = SOMADenseNDArray$new(uri, internal_use_only = "allowed_use"))
        expect_true(inherits(sdf, "SOMAArrayBase"))

        switch(tc,
               data.frame = sdf$read(iterated = TRUE),
               sparse = sdf$read_arrow_table(iterated = TRUE),
               dense = sdf$read_arrow_table(iterated = TRUE))

        expect_false(sdf$read_complete())
        dat <- sdf$read_next()
        n <- dat$num_rows
        expect_true(n > 0)

        expect_false(sdf$read_complete())
        dat <- sdf$read_next()
        n <- n + dat$num_rows
        expect_true(n > 0)

        expect_false(sdf$read_complete())
        dat <- sdf$read_next()
        n <- n + dat$num_rows
        expect_true(n > 0)

        expect_false(sdf$read_complete())
        dat <- sdf$read_next()
        n <- n + dat$num_rows
        expect_true(n > 0)

        expect_equal(n, 4848644)
        expect_true(sdf$read_complete())

        rm(sdf)
    }

})

test_that("Iterated Interface from SOMA Sparse Matrix", {
    skip_if_not_installed("pbmc3k.tiledb")      # a Suggests: pre-package 3k PBMC data

    tdir <- tempfile()
    tgzfile <- system.file("raw-data", "soco-pbmc3k.tar.gz", package="pbmc3k.tiledb")
    untar(tarfile = tgzfile, exdir = tdir)
    uri <- file.path(tdir, "soco", "pbmc3k_processed", "ms", "RNA", "X", "data")

    sdf <- SOMASparseNDArray$new(uri, internal_use_only = "allowed_use")
    expect_true(inherits(sdf, "SOMAArrayBase"))

    sdf$read_sparse_matrix(iterated = TRUE)

    nnzRows <- function(m) { sum(Matrix::rowSums(m != 0) > 0) }
    nnzTotal <- 0
    rowsTotal <- 0
    for (i in 1:4) {
        expect_false(sdf$read_complete())
        dat <- sdf$read_next()
        nnz <- Matrix::nnzero(dat)
        expect_gt(nnz, 0)
        nnzTotal <- nnzTotal + nnz
        rowsTotal <- rowsTotal + nnzRows(dat)
        # the shard dims always match the shape of the whole sparse matrix
        expect_equal(dim(dat), as.integer(sdf$shape()))
    }
    expect_true(sdf$read_complete())

    # FIXME: TileDB-SOMA issue #1111
    # expect_equal(rowsTotal, nnzRows(sdf$read_sparse_matrix()))
    # expect_equal(nnzTotal, Matrix::nnzero(sdf$read_sparse_matrix()))
    # in fact however, the test array is dense 2638x1838 with all nonzero entries.
    expect_equal(rowsTotal, 2638)
    expect_equal(nnzTotal, 4848644)
    expect_equal(nnzTotal, prod(as.integer(sdf$shape())))

    rm(sdf)

})

test_that("Iterated Interface from SOMA Dense Matrix", {
    skip_if_not_installed("pbmc3k.tiledb")      # a Suggests: pre-package 3k PBMC data

    tdir <- tempfile()
    tgzfile <- system.file("raw-data", "soco-pbmc3k.tar.gz", package="pbmc3k.tiledb")
    untar(tarfile = tgzfile, exdir = tdir)
    uri <- file.path(tdir, "soco", "pbmc3k_processed", "ms", "RNA", "X", "data")

    sdf <- SOMADenseNDArray$new(uri, internal_use_only = "allowed_use")
    expect_true(inherits(sdf, "SOMAArrayBase"))

    sdf$read_dense_matrix(iterated = TRUE)

    expect_false(sdf$read_complete())
    dat <- sdf$read_next()
    d <- dim(dat)
    expect_equal(d[2], 1838)
    n <- d[1]
    expect_true(n > 0)

    expect_false(sdf$read_complete())
    dat <- sdf$read_next()
    d <- dim(dat)
    expect_equal(d[2], 1838)
    n <- n + d[1]
    expect_true(n > 0)

    expect_false(sdf$read_complete())
    dat <- sdf$read_next()
    d <- dim(dat)
    expect_equal(d[2], 1838)
    n <- n + d[1]
    expect_true(n > 0)

    expect_false(sdf$read_complete())
    dat <- sdf$read_next()
    d <- dim(dat)
    expect_equal(d[2], 1838)
    n <- n + d[1]
    expect_true(n > 0)

    expect_equal(n, 2638)
    expect_true(sdf$read_complete())

    rm(sdf)

})
