test_that("TileDB Embedded is available", {
    # Mainly testing this function exists, without introspecting
    # overmuch on the contents.
    triple <- tiledbsoma:::tiledb_embedded_version()
    expect_true(length(triple) == 3)
})
