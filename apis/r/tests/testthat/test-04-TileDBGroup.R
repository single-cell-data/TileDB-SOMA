test_that("Non-exist", {
  uri <- file.path(withr::local_tempdir(), "new-group")
  group <- TileDBGroup$new(uri, internal_use_only = "allowed_use")

  # Should not exist on disk until created
  expect_false(dir.exists(uri))
  expect_false(group$exists())

  # Check errors on non-existent group
  expect_error(group$get("int_column"), "Item must be open for read or write.")
  expect_error(group$length(), "Item must be open for read or write.")
  expect_error(group$open(internal_use_only = "allowed_use"), "Group does not exist.")
})

test_that("Create empty", {
  uri <- file.path(withr::local_tempdir(), "new-group")
  group <- TileDBGroup$new(uri, internal_use_only = "allowed_use")

  # Create the collection on disk
  group$create(internal_use_only = "allowed_use")
  expect_error(group$create(internal_use_only = "allowed_use"), "already exists")
  expect_true(dir.exists(uri))
  expect_true(file.exists(file.path(uri, "__group")))
  expect_true(group$exists())
  fp <- file.path(uri, "__group")
  expect_match(
    get_tiledb_object_type(group$uri, group$.__enclos_env__$private$.soma_context),
    "GROUP"
  )
  group$close()
})

test_that("Accessors for empty", {
  uri <- file.path(withr::local_tempdir(), "new-group")
  group <- TileDBGroup$new(uri, internal_use_only = "allowed_use")
  group$create(internal_use_only = "allowed_use")

  group$open(mode = "READ", internal_use_only = "allowed_use")

  expect_equal(group$length(), 0)

  # Check exporters
  expect_is(group$to_list(), "list")
  expect_length(group$to_list(), 0)
  expect_is(group$to_data_frame(), "data.frame")
  expect_equal(nrow(group$to_data_frame()), 0)
  group$close()
})

test_that("Add and remove members", {
  uri <- file.path(withr::local_tempdir(), "new-group")
  group <- TileDBGroup$new(uri, internal_use_only = "allowed_use")
  group$create(internal_use_only = "allowed_use")
  group$close()

  # Create array and subgroup in isolation but do not yet add them to the group
  a1 <- TileDBArray$new(
    uri = create_empty_test_array(file.path(uri, "a1")),
    internal_use_only = "allowed_use"
  )
  g1 <- TileDBGroup$new(
    uri = tiledb::tiledb_group_create(file.path(uri, "g1")),
    internal_use_only = "allowed_use"
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

test_that("Non-relative paths", {
  uri <- file.path(withr::local_tempdir(), "new-group")
  group <- TileDBGroup$new(uri, internal_use_only = "allowed_use")
  group$create(internal_use_only = "allowed_use")

  # Error when attempting to add a relative member that's not a subpath
  g2 <- TileDBGroup$new(
    uri = file.path(withr::local_tempdir(), "not-a-subpath"),
    internal_use_only = "allowed_use"
  )
  g2$create(internal_use_only = "allowed_use")
  expect_error(
    group$set(g2, name = "g2", relative = TRUE),
    "Unable to make relative path between URIs with no common parent"
  )

  group$close()
})

test_that("Metadata", {
  uri <- file.path(withr::local_tempdir(), "group-metadata")
  group <- TileDBGroup$new(uri, internal_use_only = "allowed_use")
  expect_error(group$set_metadata(list(int_column = "float_column")), "Item must be open for write.")

  group$create(internal_use_only = "allowed_use")

  md <- list(string_column = "qux", int_column = "float_column")
  group$open("WRITE", internal_use_only = "allowed_use") # but be open for write
  group$set_metadata(md)

  # Read all metadata while the group is still open for write
  expect_equivalent(group$get_metadata("int_column"), "float_column")
  expect_equivalent(group$get_metadata("string_column"), "qux")

  readmd <- group$get_metadata()
  expect_equivalent(readmd[["string_column"]], "qux")
  expect_equivalent(readmd[["int_column"]], "float_column")
  group$close()

  # Read all metadata while the group is open for read
  group$open(mode = "READ", internal_use_only = "allowed_use")
  readmd <- group$get_metadata()
  expect_equivalent(readmd[["string_column"]], "qux")
  expect_equivalent(readmd[["int_column"]], "float_column")

  group$close()
})

# Existence proof test via cached global context
# soma_context(config = c(vfs.s3.region = "us-west-2"))
# (grp <- TileDBGroup$new(
#   uri = 's3://cellxgene-census-public-us-west-2/cell-census/2024-07-01/soma/',
#   internal_use_only = 'allowed_use'
# ))
# grp$open(mode = 'READ', internal_use_only = 'allowed_use')
# grp$names()
