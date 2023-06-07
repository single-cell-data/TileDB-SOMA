test_that("Iterated Interface from SOMAArrayReader", {
    skip_if_not_installed("pbmc3k.tiledb")      # a Suggests: pre-package 3k PBMC data
                                                # see https://ghrr.github.io/drat/
    library(arrow)
    library(bit64)
    library(tiledb)

    tdir <- tempfile()
    tgzfile <- system.file("raw-data", "soco-pbmc3k.tar.gz", package="pbmc3k.tiledb")
    untar(tarfile = tgzfile, exdir = tdir)

    uri <- file.path(tdir, "soco", "pbmc3k_processed", "ms", "RNA", "X", "data")
    expect_true(dir.exists(uri))

    ctx <- tiledb::tiledb_ctx()
    config <- tiledb::config(ctx)
    sr <- sr_setup(uri, config = as.character(config), loglevel = "warn")
    expect_true(inherits(sr, "externalptr"))

    rl <- data.frame()
    while (!tiledbsoma:::sr_complete(sr)) {
        dat <- sr_next(sr)
        D <- soma_array_to_arrow_table(dat)
        expect_true(nrow(D) > 0)
        expect_true(is_arrow_table(D))
        rl <- rbind(rl, D$to_data_frame())
    }
    expect_true(is.data.frame(rl))
    expect_equal(nrow(rl), 4848644)
    expect_equal(ncol(rl), 3)

    sr <- sr_setup(uri, config=as.character(config), dim_points=list(soma_dim_0=as.integer64(1)))
    expect_true(inherits(sr, "externalptr"))

    rl <- data.frame()
    while (!tiledbsoma:::sr_complete(sr)) {
        dat <- sr_next(sr)
        D <- soma_array_to_arrow_table(dat)
        expect_true(nrow(D) > 0)
        expect_true(is_arrow_table(D))
        rl <- rbind(rl, as.data.frame(D))
    }
    expect_true(is.data.frame(rl))
    expect_equal(nrow(rl), 1838)
    expect_equal(ncol(rl), 3)

    sr <- sr_setup(uri, config=as.character(config), dim_range=list(soma_dim_1=cbind(as.integer64(1),as.integer64(2))))
    expect_true(inherits(sr, "externalptr"))

    rl <- data.frame()
    while (!tiledbsoma:::sr_complete(sr)) {
        dat <- sr_next(sr)
        D <- soma_array_to_arrow_table(dat)
        expect_true(nrow(D) > 0)
        expect_true(is_arrow_table(D))
        rl <- rbind(rl, as.data.frame(D))
    }
    expect_true(is.data.frame(rl))
    expect_equal(nrow(rl), 5276)
    expect_equal(ncol(rl), 3)

    ## test completeness predicate on shorter data
    uri <- extract_dataset("soma-dataframe-pbmc3k-processed-obs")
    sr <- sr_setup(uri, config=as.character(config))

    expect_false(tiledbsoma:::sr_complete(sr))
    dat <- sr_next(sr)
    expect_true(tiledbsoma:::sr_complete(sr))
})


test_that("Iterated Interface from SOMA Classes", {
    skip_if_not_installed("pbmc3k.tiledb")      # a Suggests: pre-package 3k PBMC data

    tdir <- tempfile()
    tgzfile <- system.file("raw-data", "soco-pbmc3k.tar.gz", package="pbmc3k.tiledb")
    untar(tarfile = tgzfile, exdir = tdir)
    uri <- file.path(tdir, "soco", "pbmc3k_processed", "ms", "raw", "X", "data")

    ## parameterize test
    test_cases <- c("data.frame", "sparse")

    for (tc in test_cases) {
        sdf <- switch(tc,
                      data.frame = SOMADataFrame$new(uri, internal_use_only = "allowed_use"),
                      sparse = SOMASparseNDArray$new(uri, internal_use_only = "allowed_use"))
        expect_true(inherits(sdf, "SOMAArrayBase"))
        sdf$open("READ", internal_use_only = "allowed_use")

        iterator <- switch(tc,
                           data.frame = sdf$read(),
                           sparse = sdf$read()$tables())

        expect_true(inherits(iterator, "ReadIter"))
        expect_true(inherits(iterator, "TableReadIter"))

        # Test $concat()
        expect_false(iterator$read_complete())
        dat <- iterator$concat()
        expect_true(iterator$read_complete())
        expect_true(inherits(dat, "Table"))
        expect_equal(dat$num_columns, 3)
        expect_equal(dat$num_rows, 2238732)

        rm(iterator)

        # Test $read_next()
        iterator <- switch(tc,
                           data.frame = sdf$read(),
                           sparse = sdf$read()$tables())

        expect_false(iterator$read_complete())
        for (i in 1:2) {

            expect_false(iterator$read_complete())
            dat_slice <- iterator$read_next()
            expect_true(inherits(dat_slice, "Table"))
            expect_equal(dat_slice$num_columns, 3)

            if (i < 2) {
                expect_equal(dat_slice$num_rows, 2097152)
            } else {
                expect_equal(dat_slice$num_rows, 141580)
            }
        }

        expect_true(iterator$read_complete())
        expect_warning(iterator$read_next()) # returns NULL with warning
        expect_warning(iterator$read_next()) # returns NULL with warning

        rm(sdf)
    }

})

test_that("Iterated Interface from SOMA Sparse Matrix", {
    skip_if_not_installed("pbmc3k.tiledb")      # a Suggests: pre-package 3k PBMC data

    tdir <- tempfile()
    tgzfile <- system.file("raw-data", "soco-pbmc3k.tar.gz", package="pbmc3k.tiledb")
    untar(tarfile = tgzfile, exdir = tdir)
    uri <- file.path(tdir, "soco", "pbmc3k_processed", "ms", "raw", "X", "data")

    sdf <- SOMASparseNDArray$new(uri, internal_use_only = "allowed_use")
    expect_true(inherits(sdf, "SOMAArrayBase"))
    sdf$open("READ", internal_use_only = "allowed_use")

    iterator <- sdf$read()$sparse_matrix(zero_based = T)

    nnzTotal <- 0
    rowsTotal <- 0
    for (i in 1:2) {
        expect_false(iterator$read_complete())
        dat <- iterator$read_next()$get_one_based_matrix()
        nnz <- Matrix::nnzero(dat)
        expect_gt(nnz, 0)
        nnzTotal <- nnzTotal + nnz
        # the shard dims always match the shape of the whole sparse matrix
        expect_equal(dim(dat), as.integer(sdf$shape()))
    }

    expect_true(iterator$read_complete())
    expect_warning(iterator$read_next()) # returns NULL with warning
    expect_warning(iterator$read_next()) # returns NULL with warning
    expect_equal(nnzTotal, Matrix::nnzero(sdf$read()$sparse_matrix(T)$concat()$get_one_based_matrix()))
    expect_equal(nnzTotal, 2238732)

    rm(sdf)

})
