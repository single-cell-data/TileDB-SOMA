test_that("SOMAContext plumb-through", {
  skip_if(!extended_tests())

  uri <- tempfile("new-group")
  ctx <- SOMAContext$new()

  expect_s3_class(
    group <- SOMACollectionCreate(uri, context = ctx),
    "SOMACollection"
  )
  group$close()

  uri <- tempfile("new-group")
  skip_if(
    TRUE,
    message = "Disabling tests until new-style contexts are plubmed through factories"
  )
  expect_warning(
    group <- TileDBGroup$new(
      uri,
      internal_use_only = "allowed_use",
      context = ctx
    )
  )
  group$create(internal_use_only = "allowed_use")
  group$close()
})

test_that("Existence proof: create_soma_context()", {
  skip_if(
    TRUE,
    message = "Disabling tests until new-style contexts are plumbed through factories"
  )
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
      context = create_soma_context(
        config = c(vfs.s3.region = "us-west-2", vfs.s3.no_sign_request = "true")
      )
    ),
    class = 'TileDBGroup'
  )
  on.exit(grp2$close(), add = TRUE, after = FALSE)
  expect_identical(grp2$mode(), "CLOSED")
  expect_no_condition(grp2$open(
    mode = "READ",
    internal_use_only = "allowed_use"
  ))
  expect_type(grp2$names(), "character")
})

test_that("Existence proof: SOMAContext", {
  skip_if(!extended_tests() || covr_tests())
  skip_on_ci()

  uri <- "s3://cellxgene-census-public-us-west-2/cell-census/2024-07-01/soma/"
  expect_error(SOMACollectionOpen(uri))

  ctx <- SOMAContext$new(
    config = c(vfs.s3.region = "us-west-2", vfs.s3.no_sign_request = "true")
  )
  expect_s3_class(
    grp2 <- SOMACollectionOpen(uri, context = ctx),
    class = 'SOMACollection'
  )
  on.exit(grp2$close(), add = TRUE, after = FALSE)
  expect_identical(grp2$mode(), "READ")
  expect_type(grp2$names(), "character")
})
