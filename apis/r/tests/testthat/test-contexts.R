test_that("context-create", {
  skip_if(!extended_tests() || covr_tests())

  uri <- tempfile("new-group")
  ctx <- soma_context()
  group <- TileDBGroup$new(
    uri,
    internal_use_only = "allowed_use",
    soma_context = ctx
  )
  group$create(internal_use_only = "allowed_use")
  group$close()

  # Create array and subgroup in isolation but do not yet add them to the group
  a1 <- TileDBArray$new(
    uri = create_empty_test_array(file.path(uri, "a1")),
    internal_use_only = "allowed_use",
    soma_context = ctx
  )
  g1 <- TileDBGroup$new(
    uri = tiledb::tiledb_group_create(file.path(uri, "g1")),
    internal_use_only = "allowed_use",
    soma_context = ctx
  )

  # Objects are present but not yet members
  group$open(mode = "READ", internal_use_only = "allowed_use")
  expect_true(a1$exists())
  expect_true(g1$exists())
  expect_equal(group$length(), 0)
  group$close()

  # Add array and subgroup as members
  group$open(mode = "WRITE", internal_use_only = "allowed_use")
  group$set(a1, name = "a1")
  expect_equal(group$length(), 1)
  expect_equal(group$to_data_frame()$type, "ARRAY")

  group$set(g1, name = "g1")
  expect_equal(group$length(), 2)
  expect_setequal(group$to_data_frame()$type, c("ARRAY", "GROUP"))
  group$close()

  # Read back the members
  group$open(mode = "READ", internal_use_only = "allowed_use")
  expect_equal(group$length(), 2)
  expect_setequal(group$names(), c("a1", "g1"))

  # Retrieve
  o <- group$get("a1")

  expect_is(group$get("a1"), "TileDBArray")
  expect_is(group$get("g1"), "TileDBGroup")
  group$close()

  # Remove
  group$open(mode = "WRITE", internal_use_only = "allowed_use")
  group$remove("a1")
  expect_equal(group$length(), 1)
  group$remove("g1")
  expect_equal(group$length(), 0)
  group$close()

  # Remove
  group$open(mode = "READ", internal_use_only = "allowed_use")
  expect_equal(group$length(), 0)
  group$close()
})

test_that("context-fly", {
  skip_if(!extended_tests() || covr_tests())

  uri <- tempfile("new-group")
  group <- TileDBGroup$new(
    uri,
    internal_use_only = "allowed_use",
    soma_context = soma_context()
  )
  group$create(internal_use_only = "allowed_use")
  group$close()

  # Create array and subgroup in isolation but do not yet add them to the group
  a1 <- TileDBArray$new(
    uri = create_empty_test_array(file.path(uri, "a1")),
    internal_use_only = "allowed_use",
    soma_context = soma_context()
  )
  g1 <- TileDBGroup$new(
    uri = tiledb::tiledb_group_create(file.path(uri, "g1")),
    internal_use_only = "allowed_use",
    soma_context = soma_context()
  )

  # Objects are present but not yet members
  group$open(mode = "READ", internal_use_only = "allowed_use")
  expect_true(a1$exists())
  expect_true(g1$exists())
  expect_equal(group$length(), 0)
  group$close()

  # Add array and subgroup as members
  group$open(mode = "WRITE", internal_use_only = "allowed_use")
  group$set(a1, name = "a1")
  expect_equal(group$length(), 1)
  expect_equal(group$to_data_frame()$type, "ARRAY")

  group$set(g1, name = "g1")
  expect_equal(group$length(), 2)
  expect_setequal(group$to_data_frame()$type, c("ARRAY", "GROUP"))
  group$close()

  # Read back the members
  group$open(mode = "READ", internal_use_only = "allowed_use")
  expect_equal(group$length(), 2)
  expect_setequal(group$names(), c("a1", "g1"))

  # Retrieve
  o <- group$get("a1")

  expect_is(group$get("a1"), "TileDBArray")
  expect_is(group$get("g1"), "TileDBGroup")
  group$close()

  # Remove
  group$open(mode = "WRITE", internal_use_only = "allowed_use")
  group$remove("a1")
  expect_equal(group$length(), 1)
  group$remove("g1")
  expect_equal(group$length(), 0)
  group$close()

  # Remove
  group$open(mode = "READ", internal_use_only = "allowed_use")
  expect_equal(group$length(), 0)
  group$close()
})

test_that("SOMATileDBContext plumb-through", {
  skip_if(!extended_tests())

  uri <- tempfile("new-group")
  ctx <- SOMATileDBContext$new()

  expect_no_condition(group <- TileDBGroup$new(
    uri,
    internal_use_only = "allowed_use",
    tiledbsoma_ctx = ctx
  ))
  group$create(internal_use_only = "allowed_use")
  group$close()

  uri <- tempfile("new-group")
  expect_warning(group <- TileDBGroup$new(
    uri,
    internal_use_only = "allowed_use",
    tiledbsoma_ctx = ctx,
    soma_context = soma_context()
  ))
  group$create(internal_use_only = "allowed_use")
  group$close()
})

test_that("Existence proof: soma_context()", {
  skip_if(!extended_tests() || covr_tests())
  skip_on_ci()

  uri <- "s3://cellxgene-census-public-us-west-2/cell-census/2024-07-01/soma/"
  expect_s3_class(
    grp1 <- TileDBGroup$new(uri, internal_use_only = "allowed_use"),
    class = 'TileDBGroup'
  )
  on.exit(grp1$close(), add = TRUE, after = FALSE)
  expect_error(grp1$names())
  grp1$close()

  expect_s3_class(
    grp2 <- TileDBGroup$new(
      uri,
      internal_use_only = "allowed_use",
      soma_context = soma_context(config = c(vfs.s3.region = "us-west-2"))
    ),
    class = 'TileDBGroup'
  )
  on.exit(grp2$close(), add = TRUE, after = FALSE)
  expect_identical(grp2$mode(), "CLOSED")
  expect_no_condition(grp2$open(mode = "READ", internal_use_only = "allowed_use"))
  expect_type(grp2$names(), "character")
})

test_that("Existence proof: SOMATileDBContext", {
  skip_if(!extended_tests() || covr_tests())
  skip_on_ci()

  uri <- "s3://cellxgene-census-public-us-west-2/cell-census/2024-07-01/soma/"
  expect_s3_class(
    grp1 <- TileDBGroup$new(uri, internal_use_only = "allowed_use"),
    class = 'TileDBGroup'
  )
  on.exit(grp1$close(), add = TRUE, after = FALSE)
  expect_error(grp1$names())
  grp1$close()

  ctx <- SOMATileDBContext$new(config = c(vfs.s3.region = "us-west-2"))
  expect_s3_class(
    grp2 <- TileDBGroup$new(
      uri,
      internal_use_only = "allowed_use",
      tiledbsoma_ctx = ctx
    ),
    class = 'TileDBGroup'
  )
  on.exit(grp2$close(), add = TRUE, after = FALSE)
  expect_identical(grp2$mode(), "CLOSED")
  expect_no_condition(grp2$open(mode = "READ", internal_use_only = "allowed_use"))
  expect_type(grp2$names(), "character")
})
