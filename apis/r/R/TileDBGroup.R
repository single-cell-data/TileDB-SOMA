#' TileDB Group Base Class
#'
#' @description
#' Base class for interacting with TileDB groups
#' @details
#' ## Initialization
#' Upon initialization a new group is created if one does not already exist at
#' the specified `uri`.
#' @export
TileDBGroup <- R6::R6Class(
  classname = "TileDBGroup",
  inherit = TileDBObject,

  public = list(
    #' @field members Named list of members in the group
    members = list(),

    #' @description Create a new TileDBGroup object.
    #' @param uri TileDB array URI
    #' @param verbose Print status messages
    #' @param config optional configuration
    #' @param ctx optional tiledb context
    initialize = function(uri, verbose = TRUE, config = NULL, ctx = NULL) {
      super$initialize(uri, verbose, config, ctx)

      if (self$exists()) {
        if (self$verbose) {
          message(
            sprintf("Found existing %s at '%s'", self$class(), self$uri)
          )
        }
      } else {
        if (self$verbose) {
          message(
            sprintf("No %s currently exists at '%s'", self$class(), self$uri)
          )
        }
        private$create_group()
      }
      private$initialize_object()

      # Instatiate objects for existing members
      self$members <- private$instantiate_members()
    },

    #' @description Print summary of the group.
    print = function() {
      super$print()
      if (self$exists()) {
        private$format_members()
      }
    },

    #' @description Check if the group exists.
    #' @return TRUE if the group exists, FALSE otherwise.
    group_exists = function() {
      .Deprecated(
        new = "exists()",
        old = "group_exists()"
      )
      self$exists()
    },

    #' @description Return a [`tiledb_group`] object
    #' @param ... Optional arguments to pass to `tiledb::tiledb_array()`
    #' @return A [`tiledb::tiledb_group`] object.
    tiledb_group = function(...) {
      args <- list(...)
      args$uri <- self$uri
      args$ctx <- self$ctx
      do.call(tiledb::tiledb_group, args)
    },

    #' @description List the TileDB objects within the group.
    #' @param type The type of object to list, either `"ARRAY"`, or `"GROUP"`.
    #' By default all object types are listed.
    #' @return A `data.frame` with columns `URI` and `TYPE`.
    list_objects = function(type = NULL) {
      objects <- tiledb::tiledb_object_ls(self$uri)
      private$filter_by_type(objects, type)
    },

    #' @description List URIs for TileDB objects within the group.
    #' @param type The type of object to list, either `"ARRAY"`, or `"GROUP"`.
    #' By default all object types are listed.
    #' @param prefix Filter URIs whose basename contain an optional prefix.
    #' @return A character vector of object URIs with names corresponding to the
    #' basename of the object.
    list_object_uris = function(type = NULL, prefix = NULL) {
      uris <- self$list_objects(type = type)$URI
      if (is_empty(uris)) return(uris)
      names(uris) <- basename(uris)
      if (!is.null(prefix)) {
        stopifnot(is_scalar_character(prefix))
        uris <- uris[string_starts_with(names(uris), prefix)]
      }
      uris
    },

    #' @description Add new member to the group.
    #' @param object The `TileDBArray` or `TileDBGroup` object to add.
    #' @param name The name to use for the member. By default the base name of
    #' the object's URI is used.
    #' @param relative A logical value indicating whether the new member's URI
    #' is relative to the group's URI.
    #' @importFrom fs path_rel
    add_member = function(object, name = NULL, relative = NULL) {
      stopifnot(
        "Only 'TileDBArray' or 'TileDBGroup' objects can be added" =
          inherits(object, "TileDBGroup") || inherits(object, "TileDBArray"),
        is.null(name) || is_scalar_character(name)
      )

      if (is.null(relative)) {
        relative = TRUE
        if (startsWith(object$uri, "tiledb://")) {
          relative = FALSE
        }
      }

      # Because object$uri will always return an absolute URI, we need to
      # make it relative to the group's URI before adding it to the group
      if (relative) {
        uri <- fs::path_rel(object$uri, start = self$uri)
      } else {
        uri <- object$uri
      }
      name <- name %||% basename(uri)

      on.exit(private$close())
      private$open("WRITE")
      tiledb::tiledb_group_add_member(
        grp = self$object,
        uri = uri,
        relative = relative,
        name = name
      )
      self$members[[name]] <- object
    },

    #' @description Remove member from the group.
    #' @param name The name of the member to remove.
    #' @export
    remove_member = function(name) {
      on.exit(private$close())
      private$open("WRITE")
      tiledb::tiledb_group_remove_member(
        grp = self$object,
        uri = name
      )
    },

    #' @description Count the number of members in the group.
    #' @return Integer count of members in the group.
    count_members = function() {
      on.exit(private$close())
      private$open("READ")
      tiledb::tiledb_group_member_count(self$object)
    },

    #' @description List the members of the group.
    #' @param type The type of member to list, either `"ARRAY"`, or `"GROUP"`.
    #' By default all member types are listed.
    #' @return A `data.frame` with columns `URI`, `TYPE`, and `NAME`.
    list_members = function(type = NULL) {
      count <- self$count_members()
      members <- data.frame(
        TYPE = character(count),
        URI = character(count),
        NAME = character(count)
      )
      if (count == 0) return(members)

      on.exit(private$close())
      private$open("READ")
      member_list <- lapply(
        X = seq_len(count) - 1L,
        FUN = tiledb::tiledb_group_member,
        grp = self$object
      )

      members$TYPE <- vapply_char(member_list, FUN = getElement, name = 1L)
      members$URI <- vapply_char(member_list, FUN = getElement, name = 2L)
      members$NAME <- vapply_char(member_list, FUN = getElement, name = 3L)
      private$filter_by_type(members, type)
    },

    #' @description List URIs for group members
    #' @param type The type of member to list, either `"ARRAY"`, or `"GROUP"`.
    #' By default all member types are listed.
    #' @param prefix Filter for members whose name contains an optional prefix.
    #' @return A character vector of member URIs, named for the group member
    list_member_uris = function(type = NULL, prefix = NULL) {
      members <- self$list_members(type = type)
      if (is_empty(members)) return(character(0L))
      uris <- setNames(members$URI, members$NAME)
      if (!is.null(prefix)) {
        stopifnot(is_scalar_character(prefix))
        uris <- uris[string_starts_with(names(uris), prefix)]
      }
      uris
    },

    #' @description Retrieve arrays within the group that meet the specified
    #' criteria.
    #' @param type The type of group members to list, either `"ARRAY"`, or
    #' `"GROUP"`.
    #' @param prefix String prefix to filter the member names.
    #' @returns A named list of group members.
    # TODO: Add support for filtering by array metadata
    get_members = function(type = NULL, prefix = NULL) {
      if (is.null(prefix) && is.null(type)) return(self$members)
      stopifnot(
        is.null(prefix) || is_scalar_character(prefix)
      )

      matched_type <- matched_prefix <- rep(TRUE, self$count_members())
      if (!is.null(type)) {
        all_uris <- vapply_char(self$members, function(x) x$uri)
        matched_uris <- self$list_members(type = type)$URI
        matched_type <- all_uris %in% matched_uris
      }

      if (!is.null(prefix)) {
        matched_prefix <- string_starts_with(names(self$members), prefix)
      }

      self$members[matched_type & matched_prefix]
    },

    #' @description Retrieve a group member.
    #' @param name The name of the array to retrieve.
    #' @returns The array object.
    get_member = function(name) {
      stopifnot(is_scalar_character(name))
      self$members[[name]]
    },

    #' @description Retrieve metadata from the TileDB group.
    #' @param key The name of the metadata attribute to retrieve.
    #' @param prefix Filter metadata using an optional prefix. Ignored if `key`
    #'   is not NULL.
    #' @return A list of metadata values.
    get_metadata = function(key = NULL, prefix = NULL) {
      on.exit(private$close())
      private$open("READ")
      if (!is.null(key)) {
        metadata <- tiledb::tiledb_group_get_metadata(self$object, key)
      } else {
        metadata <- tiledb::tiledb_group_get_all_metadata(self$object)
        if (!is.null(prefix)) {
          metadata <- metadata[string_starts_with(names(metadata), prefix)]
        }
      }
      metadata
    },

    #' @description Add list of metadata to the TileDB group.
    #' @param metadata Named list of metadata to add.
    #' @param prefix Optional prefix to add to the metadata attribute names.
    #' @return NULL
    add_metadata = function(metadata, prefix = "") {
      stopifnot(
        "Metadata must be a named list" = is_named_list(metadata)
      )
      on.exit(private$close())
      private$open("WRITE")
      mapply(
        FUN = tiledb::tiledb_group_put_metadata,
        key = paste0(prefix, names(metadata)),
        val = metadata,
        MoreArgs = list(grp = self$object),
        SIMPLIFY = FALSE
      )
    }
  ),

  private = list(

    create_group = function() {
      if (self$verbose) {
        message(sprintf("Creating new %s at '%s'", self$class(), self$uri))
      }
      tiledb::tiledb_group_create(self$uri, ctx = self$ctx)
      private$write_object_type_metadata()
    },

    open = function(mode) {
      mode <- match.arg(mode, c("READ", "WRITE"))
      invisible(tiledb::tiledb_group_open(self$object, type = mode))
    },

    close = function() {
      invisible(tiledb::tiledb_group_close(self$object))
    },

    initialize_object = function() {
      private$tiledb_object <- tiledb::tiledb_group(self$uri, ctx = self$ctx)
      private$close()
    },

    write_object_type_metadata = function() {
      meta <- list()
      meta[[SOMA_OBJECT_TYPE_METADATA_KEY]] <- class(self)[1]
      meta[[SOMA_ENCODING_VERSION_METADATA_KEY]] <- SOMA_ENCODING_VERSION
      self$add_metadata(meta) # TileDBArray or TileDBGroup
    },

    # Instantiate existing group members
    #
    # Responsible for calling the appropriate R6 class generator for each
    # pre-existing member of a group during initialization. Currently each
    # object is named using the base name of the member's URI.
    instantiate_members = function() {
      members <- self$list_members()
      member_objects <- list()
      if (!is_empty(members)) {
        member_uris <- split(setNames(members$URI, members$NAME), members$TYPE)
        member_objects <- c(
          lapply(member_uris$ARRAY, TileDBArray$new, verbose = self$verbose),
          lapply(member_uris$GROUP, TileDBGroup$new, verbose = self$verbose)
        )
      }
      member_objects
    },

    # Filter data.frame of group objects/members by the `TYPE` column
    filter_by_type = function(x, type) {
      stopifnot(is.data.frame(x) && "TYPE" %in% names(x))
      if (is.null(type)) return(x)
      type <- match.arg(type, c("ARRAY", "GROUP"), several.ok = TRUE)
      x[x$TYPE %in% type, , drop = FALSE]
    },

    format_members = function() {
      members <- self$list_members()
      if (!is_empty(members)) {
        # denote remote URIs with *
        formatted <- paste0(
          members$NAME,
          ifelse(is_remote_uri(members$URI), "*", "")
        )
        # list by type
        formatted <- split(formatted, members$TYPE)
        if (!is.null(formatted$ARRAY)) {
          cat("  arrays:", string_collapse(sort(formatted$ARRAY)), "\n")
        }
        if (!is.null(formatted$GROUP)) {
          cat("  groups:", string_collapse(sort(formatted$GROUP)), "\n")
        }
      }
    }
  )
)
