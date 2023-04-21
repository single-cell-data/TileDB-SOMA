#' TileDB Group Base Class
#'
#' @description
#' Base class for interacting with TileDB groups (lifecycle: experimental)
#' @importFrom spdl info debug
#' @export
TileDBGroup <- R6::R6Class(
  classname = "TileDBGroup",
  inherit = TileDBObject,

  public = list(

    #' @description Print summary of the group. (lifecycle: experimental)
    print = function() {
      super$print()
      if (self$exists()) private$format_members()
    },

    #' @description Creates the data structure on disk/S3/cloud. (lifecycle: experimental)
    create = function() {
      spdl::info("Creating new {} at '{}'", self$class(), self$uri)
      tiledb::tiledb_group_create(
        uri = self$uri,
        ctx = self$get_tiledb_config('create')$context()
      )
      self
    },

    #' @description Add new member to the group. (lifecycle: experimental)
    #' @param object A `TileDBArray` or `TileDBGroup` object to add.
    #' @param name Name to use for the member. By default the base name of
    #' the object's URI is used.
    #' @param relative An optional logical value indicating whether the new
    #' object's URI is relative to the group's URI. If `NULL` (the
    #' default), the object's URI is assumed to be relative unless it is a
    #' `tiledb://` URI.
    set = function(object, name = NULL, relative = NULL) {
      stopifnot(
        "Only 'TileDBArray' or 'TileDBGroup' objects can be added" =
          inherits(object, "TileDBGroup") || inherits(object, "TileDBArray"),
        is.null(name) || is_scalar_character(name),
        is.null(relative) || is_scalar_logical(relative)
      )

      if (is.null(relative)) {
        relative <- !startsWith(object$uri, "tiledb://")
      }

      # Because object$uri will always return an absolute URI, we need to
      # make it relative to the collection's URI before adding it
      if (relative) {
        uri <- make_uri_relative(object$uri, self$uri)
      } else {
        uri <- object$uri
      }
      name <- name %||% basename(uri)

      private$open("WRITE")
      tiledb::tiledb_group_add_member(
        grp = self$object,
        uri = uri,
        relative = relative,
        name = name
      )
      # TODO: Avoid closing/re-opening the group to update the cache
      self$close()

      # We manually add the new member to member_cache in order to preserve the
      # original URI. Otherwise TileDB Cloud creation URIs are retrieved from
      # using tiledb_group_member() in the form tiledb://namespace/uuid. In this
      # form it's not possible to append new children, which is necessary during
      # ingestion.
      if (is.null(private$member_cache)) private$member_cache <- list()
      private$member_cache[[name]] <- list(
        type = tiledb::tiledb_object_type(object$uri),
        uri = object$uri,
        name = name
      )

      # We still need to update member_cache to pick-up existing members.
      # Otherwise if you open a group with existing members and add a new
      # member, the initially empty member_cache will only contain the new
      # member.
      private$update_member_cache()
    },

    #' @description Retrieve a group member by name. (lifecycle: experimental)
    #' @param name The name of the member.
    #' @returns A `TileDBArray` or `TileDBGroup`.
    get = function(name) {
      stopifnot(is_scalar_character(name))
      if (is.null(private$member_cache)) private$update_member_cache()
      member <- private$member_cache[[name]]
      if (is.null(member)) {
        stop(sprintf("No member named '%s' found", name), call. = FALSE)
      }
      private$construct_member(member$uri, member$type)
    },

    #' @description Remove member. (lifecycle: experimental)
    #' @param name Name of the member to remove.
    #' @export
    remove = function(name) {
      stopifnot(is_scalar_character(name))

      private$open("WRITE")
      on.exit(self$close())
      tiledb::tiledb_group_remove_member(
        grp = self$object,
        uri = name
      )

      # Drop member if cache has been initialized
      if (is.list(private$member_cache)) {
        private$member_cache[[name]] <- NULL
      }
    },

    #' @description Length in the number of members. (lifecycle: experimental)
    #' @return Scalar `integer`
    length = function() {
      if (is.null(private$member_cache)) private$update_member_cache()
      length(private$member_cache)
    },

    #' @description Retrieve the names of members. (lifecycle: experimental)
    #' @return A `character` vector of member names.
    names = function() {
      if (is.null(private$member_cache)) private$update_member_cache()
      names(private$member_cache) %||% character(length = 0L)
    },

    #' @description Retrieve a `list` of members. (lifecycle: experimental)
    to_list = function() {
      if (is.null(private$member_cache)) private$update_member_cache()
      private$member_cache
    },

    #' @description Retrieve a `data.frame` of members. (lifecycle: experimental)
    to_data_frame = function() {
      count <- self$length()
      df <- data.frame(
        name = character(count),
        uri = character(count),
        type = character(count)
      )
      if (count == 0) return(df)

      member_list <- self$to_list()
      df$name <- vapply_char(member_list, FUN = getElement, name = "name")
      df$uri <- vapply_char(member_list, FUN = getElement, name = "uri")
      df$type <- vapply_char(member_list, FUN = getElement, name = "type")
      df
    },

    #' @description Retrieve metadata. (lifecycle: experimental)
    #' @param key The name of the metadata attribute to retrieve.
    #'   is not NULL.
    #' @return A list of metadata values.
    get_metadata = function(key = NULL) {
      on.exit(self$close())
      private$open("READ")
      spdl::debug("Retrieving metadata for {} '{}'", self$class(), self$uri)
      if (!is.null(key)) {
        return(tiledb::tiledb_group_get_metadata(self$object, key))
      } else {
        return(tiledb::tiledb_group_get_all_metadata(self$object))
      }
    },

    #' @description Add list of metadata. (lifecycle: experimental)
    #' @param metadata Named list of metadata to add.
    #' @return NULL
    set_metadata = function(metadata) {
      stopifnot(
        "Metadata must be a named list" = is_named_list(metadata)
      )
      private$open("WRITE")
      on.exit(self$close())
      spdl::debug("Writing metadata to {} '{}'", self$class(), self$uri)
      dev_null <- mapply(
        FUN = tiledb::tiledb_group_put_metadata,
        key = names(metadata),
        val = metadata,
        MoreArgs = list(grp = self$object),
        SIMPLIFY = FALSE
      )
    },

    #' @description Close the SOMA object.
    #' @return The object, invisibly
    close = function() {
      spdl::debug("Closing {} '{}'", self$class(), self$uri)
      invisible(tiledb::tiledb_group_close(self$object))
    }

  ),

  private = list(

    # @description List of cached group members
    # Initially NULL, once the group is created or opened, this is populated
    # with a list that's empty or contains the group members.
    member_cache = NULL,

    open = function(mode) {
      mode <- match.arg(mode, c("READ", "WRITE"))
      spdl::debug(
        "Opening {} '{}' in {} mode", self$class(), self$uri, mode
      )
      invisible(tiledb::tiledb_group_open(self$object, type = mode))
    },

    initialize_object = function() {
      private$tiledb_object <- tiledb::tiledb_group(
        self$uri,
        ctx = self$tiledbsoma_ctx$context()
      )
      self$close()
    },

    # @description Retrieve all group members. (lifecycle: experimental)
    # @return A list indexed by group member names where each element is a
    # list with names: name, uri, and type.
    get_all_members = function() {
      private$open("READ")
      on.exit(self$close())

      count <- tiledb::tiledb_group_member_count(self$object)
      if (count == 0) return(list())

      members <- vector(mode = "list", length = count)
      if (count == 0) return(members)

      for (i in seq_len(count)) {
        members[[i]] <- setNames(
          object = as.list(tiledb::tiledb_group_member(self$object, i - 1L)),
          nm = c("type", "uri", "name")
        )
      }

      # Key the list by group member name
      names(members) <- vapply_char(members, FUN = getElement, name = "name")
      members
    },

    update_member_cache = function() {
      spdl::debug("Updating member cache for {} '{}'", self$class(), self$uri)
      members <- private$get_all_members()
      if (is.null(private$member_cache)) {
        private$member_cache <- members
      } else {
        # Don't clobber existing cache members in order to retain original URIs
        members <- members[setdiff(names(members), names(private$member_cache))]
        private$member_cache <- utils::modifyList(private$member_cache, members)
      }
    },

    # Instantiate a group member object.
    # Responsible for calling the appropriate R6 class constructor.
    construct_member = function(uri, type) {
      stopifnot(
        is_scalar_character(uri),
        is_scalar_character(type)
      )
      constructor <- switch(type,
        ARRAY = TileDBArray$new,
        GROUP = TileDBGroup$new,
        stop(sprintf("Unknown member type: %s", type), call. = FALSE)
      )
      constructor(uri, tiledbsoma_ctx = self$tiledbsoma_ctx,
                  platform_config = self$platform_config, internal_use_only = "allowed_use")
    },

    format_members = function() {
      members <- self$to_data_frame()

      if (nrow(members) > 0) {
        # denote remote URIs with *
        formatted <- paste0(
          members$name,
          ifelse(is_remote_uri(members$uri), "*", "")
        )
        # list by type
        formatted <- split(formatted, members$type)
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
