test_that("Iterated Interface from SOMAReader", {
    skip_if_not_installed("pbmc3k.tiledb")      # a Suggests: pre-package 3k PBMC data

    library(arch)
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
    sr <- sr_setup(ctx@ptr, uri)
    expect_true(inherits(sr, "externalptr"))
    rl <- data.frame()
    while (nrow(rl) == 0 || !tiledbsoma:::sr_complete(sr)) {
        dat <- sr_next(sr)
        D <- tiledbsoma:::arrow_to_dt(dat)
        expect_true(nrow(D) > 0)
        expect_true(inherits(D, "data.table"))
        rl <- rbind(rl, D)
    }
    expect_true(inherits(rl, "data.table"))
    expect_equal(nrow(rl), 4848644)
    expect_equal(ncol(rl), 3)

    sr <- sr_setup(ctx@ptr, uri, dim_points=list(soma_dim_0=as.integer64(1)))
    expect_true(inherits(sr, "externalptr"))
    rl <- data.frame()
    while (nrow(rl) == 0 || !tiledbsoma:::sr_complete(sr)) {
        dat <- sr_next(sr)
        D <- tiledbsoma:::arrow_to_dt(dat)
        expect_true(nrow(D) > 0)
        expect_true(inherits(D, "data.table"))
        rl <- rbind(rl, D)
    }
    expect_true(inherits(rl, "data.table"))
    expect_equal(nrow(rl), 1838)
    expect_equal(ncol(rl), 3)

    sr <- sr_setup(ctx@ptr, uri, dim_range=list(soma_dim_1=cbind(as.integer64(1),as.integer64(2))))
    expect_true(inherits(sr, "externalptr"))
    rl <- data.frame()
    while (nrow(rl) == 0 || !tiledbsoma:::sr_complete(sr)) {
        dat <- sr_next(sr)
        D <- tiledbsoma:::arrow_to_dt(dat)
        expect_true(nrow(D) > 0)
        expect_true(inherits(D, "data.table"))
        rl <- rbind(rl, D)
    }
    expect_true(inherits(rl, "data.table"))
    expect_equal(nrow(rl), 5276)
    expect_equal(ncol(rl), 3)
})
