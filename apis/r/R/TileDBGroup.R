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

    mode = function() {
      private$.mode
    },

    #' @description Open the SOMA object for read or write.
    #' @param internal_use_only Character value to signal 'permitted' call as
    #' `new()` is considered internal and should not be called directly
    #' @return The object, invisibly
    open = function(mode="READ", internal_use_only = NULL) {
      mode <- match.arg(mode, c("READ", "WRITE"))
      if (is.null(internal_use_only) || internal_use_only != "allowed_use") {
        stop(paste("Use of the open() method is discouraged. Consider using a",
                   "factory method as e.g. 'SOMADataFrameOpen()'."), call. = FALSE)
      }
      spdl::debug(
        "Opening {} '{}' in {} mode", self$class(), self$uri, mode
      )
      private$.mode = mode

      private$.tiledb_group <- tiledb::tiledb_group(
        self$uri,
        type = mode,
        ctx = private$.tiledb_ctx
      )

      private$update_member_cache()
      # XXX TODO
      #private$update_metadata_cache()

      invisible(private$.tiledb_group)
    },

    #' @description Close the SOMA object.
    #' @return The object, invisibly
    close = function() {
      private$check_open_for_read_or_write()

      spdl::debug("Closing {} '{}'", self$class(), self$uri)
      invisible(tiledb::tiledb_group_close(private$.tiledb_group))
      private$.mode <- "CLOSED"
      private$.tiledb_group <- NULL
    },

    #' @description Creates the data structure on disk/S3/cloud. (lifecycle: experimental)
    #' @param internal_use_only Character value to signal 'permitted' call as
    #' `new()` is considered internal and should not be called directly
    create = function(internal_use_only = NULL) {
      if (is.null(internal_use_only) || internal_use_only != "allowed_use") {
        stop(paste("Use of the create() method is discouraged. Consider using a",
                   "factory method as e.g. 'SOMACollectionCreate()'."), call. = FALSE)
      }

      spdl::info("Creating new {} at '{}'", self$class(), self$uri)
      tiledb::tiledb_group_create(uri = self$uri, ctx = private$.tiledb_ctx)

      self$open("WRITE", internal_use_only = "allowed_use")

      invisible(private$.tiledb_group)
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

      private$check_open_for_write()

      tiledb::tiledb_group_add_member(
        grp = private$.tiledb_group,
        uri = uri,
        relative = relative,
        name = name
      )

      # XXX EXTRACT METHOD

      # We explicitly add the new member to member_cache in order to preserve the
      # original URI. Otherwise TileDB Cloud creation URIs are retrieved from
      # using tiledb_group_member() in the form tiledb://namespace/uuid. In this
      # form it's not possible to append new children, which is necessary during
      # ingestion.
      if (is.null(private$.member_cache)) private$.member_cache <- list()
      private$.member_cache[[name]] <- list(
        type = tiledb::tiledb_object_type(object$uri),
        uri = object$uri,
        name = name
      )

      # We still need to update member_cache to pick up existing members.
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
      private$check_open_for_read_or_write()

      private$fill_member_cache_if_null()
      member <- private$.member_cache[[name]]
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
      private$check_open_for_write()

      tiledb::tiledb_group_remove_member(
        grp = private$.tiledb_group,
        uri = name
      )

      # Drop member if cache has been initialized
      if (is.list(private$.member_cache)) {
        private$.member_cache[[name]] <- NULL
      }
    },

    #' @description Length in the number of members. (lifecycle: experimental)
    #' @return Scalar `integer`
    length = function() {
      private$check_open_for_read_or_write()
      private$fill_member_cache_if_null()
      length(private$.member_cache)
    },

    #' @description Retrieve the names of members. (lifecycle: experimental)
    #' @return A `character` vector of member names.
    names = function() {
      private$check_open_for_read_or_write()
      private$fill_member_cache_if_null()
      names(private$.member_cache) %||% character(length = 0L)
    },

    #' @description Retrieve a `list` of members. (lifecycle: experimental)
    to_list = function() {
      private$check_open_for_read_or_write()
      private$fill_member_cache_if_null()
      private$.member_cache
    },

    #' @description Retrieve a `data.frame` of members. (lifecycle: experimental)
    to_data_frame = function() {
      private$check_open_for_read_or_write()
      count <- self$length()
      df <- data.frame(
        name = character(count),
        uri  = character(count),
        type = character(count)
      )
      if (count == 0) return(df)

      member_list <- self$to_list()
      df$name <- vapply_char(member_list, FUN = getElement, name = "name")
      df$uri <- vapply_char(member_list, FUN = getElement, name = "uri")
      df$type <- vapply_char(member_list, FUN = getElement, name = "type")
      df
    },

    # XXX MAKE A CACHING LAYER FOR METADATA TOO

    #' @description Retrieve metadata. (lifecycle: experimental)
    #' @param key The name of the metadata attribute to retrieve.
    #'   is not NULL.
    #' @return A list of metadata values.
    get_metadata = function(key = NULL) {
      private$check_open_for_read_or_write()

      spdl::debug("Retrieving metadata for {} '{}'", self$class(), self$uri)
      if (!is.null(key)) {
        return(tiledb::tiledb_group_get_metadata(private$.tiledb_group, key))
      } else {
        return(tiledb::tiledb_group_get_all_metadata(private$.tiledb_group))
      }
    },

    #' @description Add list of metadata. (lifecycle: experimental)
    #' @param metadata Named list of metadata to add.
    #' @return NULL
    set_metadata = function(metadata) {
      stopifnot(
        "Metadata must be a named list" = is_named_list(metadata)
      )

      private$check_open_for_write()

      spdl::debug("Writing metadata to {} '{}'", self$class(), self$uri)
      dev_null <- mapply(
        FUN = tiledb::tiledb_group_put_metadata,
        key = names(metadata),
        val = metadata,
        MoreArgs = list(grp = private$.tiledb_group),
        SIMPLIFY = FALSE
      )
    }

  ),

  private = list(

    # Pro tip: in R6 we can't set these to anything other than NULL here, even if we want to.  If
    # you want them defaulted to anything other than NULL, leave them NULL here and set the defaults
    # in the constructor.
    .mode = NULL,

    # @description This is a handle at the TileDB-R level
    .tiledb_group = NULL,

    # @description List of cached group members
    # Initially NULL, once the group is created or opened, this is populated
    # with a list that's empty or contains the group members.
    .member_cache = NULL,

    # XXX TODO
    #.metadta_cache = NULL,

    check_ever_opened = function() {
      stopifnot(
        "Group has not been opened." = !is.null(private$.tiledb_group)
      )
    },

    # Per the spec, invoking user-level read requires open for read mode.
    check_open_for_read = function() {
      private$check_ever_opened()
      stopifnot(
        "Group must be open for read." = private$.mode == "READ"
      )
    },

    # Per the spec, invoking user-level write requires open for read mode.
    check_open_for_write = function() {
      private$check_ever_opened()
      stopifnot(
        "Group must be open for write." = private$.mode == "WRITE"
      )
    },

    # Per the spec, invoking user-level get-metadata requires open for read mode or write mode.
    check_open_for_read_or_write = function() {
      private$check_ever_opened()
      stopifnot(
        "Group must be open for read or write." = (private$.mode == "READ" || private$.mode == "WRITE")
      )
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

    # @description Retrieve all group members. (lifecycle: experimental)
    # @return A list indexed by group member names where each element is a
    # list with names: name, uri, and type.
    get_all_members_uncached_read = function(group_handle) {

      count <- tiledb::tiledb_group_member_count(group_handle)
      if (count == 0) return(list())

      members <- vector(mode = "list", length = count)
      if (count == 0) return(members)

      for (i in seq_len(count)) {
        members[[i]] <- setNames(
          object = as.list(tiledb::tiledb_group_member(group_handle, i - 1L)),
          nm = c("type", "uri", "name")
        )
      }

      # Key the list by group member name
      names(members) <- vapply_char(members, FUN = getElement, name = "name")
      members
    },

    fill_member_cache_if_null = function() {
      if (is.null(private$.member_cache)) {
        private$update_member_cache()
      }
    },

    update_member_cache = function() {
      spdl::debug("Updating member cache for {} '{}'", self$class(), self$uri)

      # XXX COMMENT WHY -- PASTE IN NOTES
      group_handle <- private$.tiledb_group
      if (private$.mode == "WRITE") {
        group_handle <- tiledb::tiledb_group(self$uri, type = "READ", ctx = private$.tiledb_ctx)
      }

      members <- private$get_all_members_uncached_read(group_handle)
      if (is.null(private$.member_cache)) {
        private$.member_cache <- members
      } else {
        # Don't clobber existing cache members, in order to retain original URIs
        members <- members[setdiff(names(members), names(private$.member_cache))]
        private$.member_cache <- utils::modifyList(private$.member_cache, members)
      }

      if (private$.mode == "WRITE") {
        tiledb::tiledb_group_close(group_handle)
      }
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
