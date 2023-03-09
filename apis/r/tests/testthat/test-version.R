test_that("TileDB Embedded is available", {
    # Mainly testing this function exists, without introspecting
    # overmuch on the contents.
    version_string <- tiledbsoma:::tiledb_embedded_version()
    expect_true(nchar(version_string) > 10)
})
