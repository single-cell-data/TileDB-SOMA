test_that("Basic SOMAReader", {
    library(arch)

    tdir <- tempfile()
    tgzfile <- system.file("raw-data", "soco-pbmc3k_processed-obs.tar.gz", package="tiledbsoma")
    untar(tarfile = tgzfile, exdir = tdir)

    uri <- file.path(tdir, "obs")
    names <- get_column_names(uri)
    expect_equal(length(names), 6L)
    expect_equal(names, c("soma_rowid", "obs_id", "n_genes", "percent_mito", "n_counts", "louvain"))

    for (n in names) {
        col <- get_column(uri, n)
        expect_equal(length(col), 2638L)
    }

    df <- get_table(uri)
    expect_equal(nrow(df), 2638L)
    expect_equal(ncol(df), 6L)

    tiledbsoma:::set_log_level("warn") # suppress info level
    dfsub <- export_column_direct(uri, c("n_genes", "percent_mito"))
    expect_equal(df$`n_genes`, from_arch_array(dfsub$`n_genes`))
    expect_equal(df$`percent_mito`, from_arch_array(dfsub$`percent_mito`))

    columns <- c("n_counts", "n_genes", "louvain")
    z <- export_arrow_array(uri, columns)
    expect_true(inherits(z, "arch_array"))
})
