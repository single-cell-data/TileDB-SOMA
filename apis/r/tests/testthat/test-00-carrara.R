test_that("is_relative_uri checks", {
    expect_true(is_relative_uri("foo"))
    expect_true(is_relative_uri("foo/bar"))
    expect_false(is_relative_uri("/foo/bar"))
    expect_false(is_relative_uri("file://foo/bar"))
})

test_that("sanitize_key checks", {
    # v3
    expect_equal(sanitize_key("foo", "tiledbv3"), "foo")
    # Test that it errors on slash
    expect_error(sanitize_key("foo/bar", "tiledbv3"), "must not contain slash")

    # v2
    expect_equal(sanitize_key("foo", "tiledbv2"), "foo")
    # Encoded
    expect_true(sanitize_key("foo/bar", "tiledbv2") != "foo/bar")
})

test_that("SOMATileDBContext detects data protocol", {
    skip_if_not_installed("tiledb")

    # Initialize with simple config
    ctx <- SOMATileDBContext$new()

    expect_equal(ctx$data_protocol("file://foo/bar"), "tiledbv2")
    expect_equal(ctx$data_protocol("s3://foo/bar"), "tiledbv2")
    expect_equal(ctx$data_protocol("tiledb://foo/bar"), "tiledbv3")

    # Cloud check
    ctx_cloud <- SOMATileDBContext$new(c("rest.server_address" = "https://api.tiledb.com"))
    expect_equal(ctx_cloud$data_protocol("tiledb://foo/bar"), "tiledbv2")

    # Carrara check
    ctx_carrara <- SOMATileDBContext$new(c("rest.server_address" = "https://api.something-else.com"))
    expect_equal(ctx_carrara$data_protocol("tiledb://foo/bar"), "tiledbv3")
})
