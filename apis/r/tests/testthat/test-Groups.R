test_that("Group", {
    #uri <- tempfile()
    uri <- "/tmp/tiledb/grouptest"; if (dir.exists(uri)) unlink(uri, recursive=TRUE)

    ts1 <- as.POSIXct(1, tz="UTC")
    grp <- TileDBGroup$new(uri, NULL, NULL, ts1, "allowed_use")

    ## Should not exist on disk until created
    expect_false(dir.exists(uri))
    expect_false(grp$exists())

    grp$create(internal_use_only = "allowed_use")
    expect_true(dir.exists(uri))
    expect_true(grp$exists())

    res <- grp$open("READ", internal_use_only = "allowed_use")
    #print(class(res))

    expect_equal(grp$length(), 0)

    ## Check exporters
    expect_is(grp$to_list(), "list")
    expect_length(grp$to_list(), 0)
    expect_is(grp$to_data_frame(), "data.frame")
    expect_equal(nrow(grp$to_data_frame()), 0)
    grp$close()

})


test_that("Add and remove members", {
  #uri <- file.path(withr::local_tempdir(), "new-group")
  uri <- "/tmp/tiledb/grouptest"; if (dir.exists(uri)) unlink(uri, recursive=TRUE)
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

  # Remove: check with rebuilt cache
  group$open(mode = "READ", internal_use_only = "allowed_use")
  expect_equal(group$length(), 0)
  group$close()

  #rm(group); gc()  # forces finalizer and then destructor

})

test_that("Non-relative paths", {
  #uri <- file.path(withr::local_tempdir(), "new-group")
  uri <- "/tmp/tiledb/grouptest"; if (dir.exists(uri)) unlink(uri, recursive=TRUE)
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
  #uri <- file.path(withr::local_tempdir(), "group-metadata")
  uri <- "/tmp/tiledb/grouptest"; if (dir.exists(uri)) unlink(uri, recursive=TRUE)
  group <- TileDBGroup$new(uri, internal_use_only = "allowed_use")
  expect_error(group$set_metadata(list(foo = "bar")), "Item must be open for write.")

  group$create(internal_use_only = "allowed_use")

  md <- list(baz = "qux", foo = "bar")
  group$open("WRITE", internal_use_only = "allowed_use") ## FIXME
  group$set_metadata(md)

  # Read all metadata while the group is still open for write
  expect_equivalent(group$get_metadata("foo"), "bar")
  expect_equivalent(group$get_metadata("baz"), "qux")

  readmd <- group$get_metadata()
  expect_equivalent(readmd[["baz"]], "qux")
  expect_equivalent(readmd[["foo"]], "bar")
  group$close()

  # Read all metadata while the group is open for read
  group$open(mode = "READ", internal_use_only = "allowed_use")
  readmd <- group$get_metadata()
  expect_equivalent(readmd[["baz"]], "qux")
  expect_equivalent(readmd[["foo"]], "bar")

  group$close()

})
