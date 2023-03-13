test_that("libtiledbsoma version info is available", {
    # Mainly testing this function exists, without introspecting
    # overmuch on the contents.
    triple <- tiledbsoma:::tiledb_embedded_version()
    expect_true(length(triple) == 3)

    v <- tiledbsoma:::libtiledbsoma_version()
    expect_type(v, "character")
    expect_gt(length(v), 0)
})
