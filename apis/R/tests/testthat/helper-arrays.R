create_empty_test_array <- function(uri) {
  stopifnot(!dir.exists(uri))
  dim <- tiledb::tiledb_dim("d0", type = "ASCII", domain = NULL, tile = NULL)
  dom <- tiledb::tiledb_domain(dims = dim)
  schema <- tiledb::tiledb_array_schema(
    domain = dom,
    attrs = c(tiledb::tiledb_attr("a", type = "INT32")),
    sparse = TRUE
  )
  tiledb::tiledb_array_create(uri, schema)
  return(uri)
}

# Create a parent group at the specified URI that contains a configurable number
# of sub-arrays (all of which are empty) and sub-groups.
# param n_arrays: number of arrays to create, named a1..an
# param n_groups: number of groups to create, named g1..gn
create_test_group_with_members <- function(
  uri,
  relative = TRUE,
  n_arrays = 1,
  n_groups = 1
) {
  stopifnot(
    is.logical(relative),
    n_arrays >= 0,
    n_groups >= 0
  )

  # create parent group
  tiledb::tiledb_group_create(uri)

  # create sub-arrays/-groups
  members <- character(length = 0L)

  if (n_arrays > 0) {
    array_uris <- file.path(uri, paste0("a", seq_len(n_arrays)))
    arrays <- vapply_char(
      setNames(array_uris, basename(array_uris)),
      create_empty_test_array
    )
    members <- c(members, arrays)
  }

  if (n_groups > 0) {
    group_uris <- file.path(uri, paste0("g", seq_len(n_groups)))
    groups <- vapply_char(
      setNames(group_uris, basename(group_uris)),
      tiledb::tiledb_group_create
    )
    members <- c(members, groups)
  }

  # add members to parent group
  grp <- tiledb::tiledb_group(uri, "WRITE")

  Map(function(member_uri, member_name) {
      if (relative) member_uri <- basename(member_uri)
      tiledb::tiledb_group_add_member(
        grp = grp,
        uri = member_uri,
        relative = relative,
        name = member_name
      )
    },
    member_uri = members,
    member_name = names(members)
  )

  tiledb::tiledb_group_close(grp)
  return(uri)
}
