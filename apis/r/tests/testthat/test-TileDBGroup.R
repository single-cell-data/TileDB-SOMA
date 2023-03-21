
test_that("Basic mechanics", {
  uri <- file.path(withr::local_tempdir(), "new-group")
  group <- TileDBGroup$new(uri, internal_use_only = "allowed_use")

  # Should not exist on disk until created
  expect_false(dir.exists(uri))
  expect_false(group$exists())

  # Check errors on non-existent group
  expect_error(group$get("foo"), "Group does not exist.")
  expect_error(group$length(), "Group does not exist.")

  # Create the collection on disk
  group$create()
  expect_error(group$create(), "already exists")
  expect_true(dir.exists(uri))
  expect_true(group$exists())
  expect_match(tiledb::tiledb_object_type(uri), "GROUP")
  expect_equal(group$length(), 0)

  # Check exporters
  expect_is(group$to_list(), "list")
  expect_length(group$to_list(), 0)
  expect_is(group$to_data_frame(), "data.frame")
  expect_equal(nrow(group$to_data_frame()), 0)

  # Add members to the group
  a1 <- TileDBArray$new(
    uri = create_empty_test_array(file.path(uri, "a1")),
    internal_use_only = "allowed_use"
  )
  g1 <- TileDBGroup$new(
    uri = tiledb::tiledb_group_create(file.path(uri, "g1")),
    internal_use_only = "allowed_use"
  )

  # Objects are present but not yet members
  expect_true(a1$exists())
  expect_true(g1$exists())
  expect_equal(group$length(), 0)

  # Add sub-array/group as members
  group$set(a1, name = "a1")
  expect_equal(group$length(), 1)
  expect_equal(group$to_data_frame()$type, "ARRAY")

  group$set(g1, name = "g1")
  expect_equal(group$length(), 2)
  expect_setequal(group$to_data_frame()$type, c("ARRAY", "GROUP"))

  # Read back the members
  group_readback <- TileDBGroup$new(group$uri, internal_use_only = "allowed_use")
  expect_equal(group_readback$length(), 2)
  expect_setequal(group$names(), c("a1", "g1"))

  # Retrieve
  expect_is(group_readback$get("a1"), "TileDBArray")
  expect_is(group_readback$get("g1"), "TileDBGroup")

  # Error when attempting to add a relative member that's not a subpath
  g2 <- TileDBGroup$new(
    uri = file.path(withr::local_tempdir(), "not-a-subpath"),
    internal_use_only = "allowed_use"
  )$create()
  expect_error(
    group$set(g2, name = "g2", relative = TRUE),
    "Unable to make relative path between URIs with no common parent"
  )

  # Remove
  group_readback$remove("a1")
  expect_equal(group_readback$length(), 1)
  group_readback$remove("g1")
  expect_equal(group_readback$length(), 0)
})

test_that("Metadata", {
  uri <- file.path(withr::local_tempdir(), "group-metadata")
  group <- TileDBGroup$new(uri, internal_use_only = "allowed_use")
  expect_error(group$set_metadata(list(foo = "bar")), "Group does not exist.")

  group$create()
  md <- list(baz = "qux", foo = "bar")
  group$set_metadata(md)
  expect_equivalent(group$get_metadata("foo"), "bar")
  expect_equivalent(group$get_metadata("baz"), "qux")

  # Read all metadata
  readmd <- group$get_metadata()
  expect_equivalent(readmd[["baz"]], "qux")
  expect_equivalent(readmd[["foo"]], "bar")
})
