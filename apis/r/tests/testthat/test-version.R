test_that("version info for libtiledbsoma and its linked TileDB Embedded", {
    # Mainly testing these functions exist, without introspecting
    # overmuch on the contents.
    triple <- tiledbsoma:::tiledb_embedded_version()
    expect_true(length(triple) == 3)

    v <- tiledbsoma:::libtiledbsoma_version()
    expect_type(v, "character")
    expect_gt(length(v), 0)
})
