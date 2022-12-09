#' TileDB Group Base Class
#'
#' @description
#' Base class for interacting with TileDB groups
#' @export
TileDBGroup <- R6::R6Class(
  classname = "TileDBGroup",
  inherit = TileDBObject,

  public = list(

    #' @description Print summary of the group.
    print = function() {
      super$print()
      if (self$exists()) private$format_members()
    },

    #' @description Creates the data structure on disk/S3/cloud.
    create = function() {
      message(sprintf("Creating new %s at '%s'", self$class(), self$uri))
      tiledb::tiledb_group_create(self$uri, ctx = self$ctx)
      private$write_object_type_metadata()
    },

    #' @description Add new member to the group.
    #' @param object A `TileDBArray` or `TileDBGroup` object to add.
    #' @param name Name to use for the member. By default the base name of
    #' the object's URI is used.
    #' @param relative An optional logical value indicating whether the new
    #' object's URI is relative to the group's URI. If `NULL` (the
    #' default), the object's URI is assumed to be relative unless it is a
    #' `tiledb://` URI.
    #' @importFrom fs path_rel
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
        uri <- fs::path_rel(object$uri, start = self$uri)
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
      private$close()
      private$update_member_cache()
    },

    #' @description Retrieve a group member by name.
    #' @param name The name of the member.
    #' @returns A `TileDBArray` or `TileDBGroup`.
    get = function(name) {
      stopifnot(is_scalar_character(name))
      if (is_empty(private$member_cache)) private$update_member_cache()
      member <- private$member_cache[[name]]
      if (is.null(member)) {
        stop(sprintf("No member named '%s' found", name))
      }
      private$construct_member(member$uri, member$type)
    },

    #' @description Remove member.
    #' @param name Name of the member to remove.
    #' @export
    remove = function(name) {
      stopifnot(is_scalar_character(name))

      private$open("WRITE")
      tiledb::tiledb_group_remove_member(
        grp = self$object,
        uri = name
      )
      private$close()
      # TODO: Avoid closing/re-opening the group to update the cache
      private$update_member_cache()
    },

    #' @description Length in the number of members.
    #' @return Scalar `integer`
    length = function() {
      if (is_empty(private$member_cache)) private$update_member_cache()
      length(private$member_cache)
    },

    #' @description Retrieve a `list` of members.
    to_list = function() {
      if (is_empty(private$member_cache)) private$update_member_cache()
      private$member_cache
    },

    #' @description Retrieve a `data.frame` of members.
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

    #' @description Retrieve metadata.
    #' @param key The name of the metadata attribute to retrieve.
    #'   is not NULL.
    #' @return A list of metadata values.
    get_metadata = function(key = NULL) {
      on.exit(private$close())
      private$open("READ")
      if (!is.null(key)) {
        return(tiledb::tiledb_group_get_metadata(self$object, key))
      } else {
        return(tiledb::tiledb_group_get_all_metadata(self$object))
      }
    },

    #' @description Add list of metadata.
    #' @param metadata Named list of metadata to add.
    #' @return NULL
    set_metadata = function(metadata) {
      stopifnot(
        "Metadata must be a named list" = is_named_list(metadata)
      )
      on.exit(private$close())
      private$open("WRITE")
      dev_null <- mapply(
        FUN = tiledb::tiledb_group_put_metadata,
        key = names(metadata),
        val = metadata,
        MoreArgs = list(grp = self$object),
        SIMPLIFY = FALSE
      )
    }
  ),

  private = list(

    # @description List of cached group members
    member_cache = list(),

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
      meta[[SOMA_OBJECT_TYPE_METADATA_KEY]] <- self$class()
      meta[[SOMA_ENCODING_VERSION_METADATA_KEY]] <- SOMA_ENCODING_VERSION
      self$set_metadata(meta)
    },

    # @description Retrieve all group members.
    # @return A list indexed by group member names where each element is a
    # list with names: name, uri, and type.
    get_all_members = function() {
      on.exit(private$close())

      private$open("READ")
      count <- tiledb::tiledb_group_member_count(self$object)
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
      message("Updating member cache")
      private$member_cache <- private$get_all_members()
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
        stop(sprintf("Unknown member type: %s", type))
      )
      constructor(uri, ctx = self$ctx, platform_config = self$platform_config)
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
