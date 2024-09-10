#' TileDB Group Base Class
#'
#' @description
#' Base class for interacting with TileDB groups (lifecycle: maturing)
#' @importFrom spdl info debug
#' @keywords internal
#' @export
TileDBGroup <- R6::R6Class(
  classname = "TileDBGroup",
  inherit = TileDBObject,

  public = list(

    #' @description Print summary of the group. (lifecycle: maturing)
    print = function() {
      super$print()
      if (self$exists()) {
        if (self$is_open()) {
          private$format_members()
        } else {
          cat("  closed\n")
        }
      }
    },

    #' @description Creates the data structure on disk/S3/cloud. (lifecycle: maturing)
    #' @param internal_use_only Character value to signal this is a 'permitted' call,
    #' as `create()` is considered internal and should not be called directly.
    create = function(internal_use_only = NULL) {
      if (is.null(internal_use_only) || internal_use_only != "allowed_use") {
        stop(paste("Use of the create() method is for internal use only. Consider using a",
                   "factory method as e.g. 'SOMACollectionCreate()'."), call. = FALSE)
      }

      spdl::debug("[TileDBGroup$create] Creating new {} at '{}' at {}",
                  self$class(), self$uri, self$tiledb_timestamp)

      private$.soma_context <- soma_context()  # FIXME via factory and paramater_config
      c_group_create(self$uri, self$class(), private$.soma_context,
                     self$.tiledb_timestamp_range)  ## FIXME: use to be added accessor

      invisible(private$.tiledb_group)
    },

    #' @description Open the SOMA object for read or write.
    #' @param mode Mode to open in; defaults to `READ`.
    #' @param internal_use_only Character value to signal this is a 'permitted' call,
    #' as `open()` is considered internal and should not be called directly.
    #' @return The object, invisibly
    open = function(mode=c("READ", "WRITE"), internal_use_only = NULL) {
      mode <- match.arg(mode)
      if (is.null(internal_use_only) || internal_use_only != "allowed_use") {
        stop(paste("Use of the open() method is for internal use only. Consider using a",
                   "factory method as e.g. 'SOMACollectionOpen()'."), call. = FALSE)
      }
      private$.mode = mode
      private$.group_open_timestamp <- if (mode == "READ" && is.null(self$tiledb_timestamp)) {
        # In READ mode, if the opener supplied no timestamp then we default to the time of
        # opening, providing a temporal snapshot of all group members.
        Sys.time()
      } else {
        self$tiledb_timestamp
      }
      if (is.null(private$.group_open_timestamp)) {
        spdl::debug("[TileDBGroup$open] Opening {} '{}' in {} mode", self$class(), self$uri, mode)
        private$.tiledb_group <- c_group_open(self$uri, type = mode, ctx = soma_context())#private$.some_context)
      } else {
        if (internal_use_only != "allowed_use") stopifnot("tiledb_timestamp not yet supported for WRITE mode" = mode == "READ")
        spdl::debug("[TileDBGroup$open] Opening {} '{}' in {} mode at {} ptr null {}",
                    self$class(), self$uri, mode, private$.group_open_timestamp,
                    is.null(private$.soma_context))
        ## The Group API does not expose a timestamp setter so we have to go via the config
        private$.tiledb_group <- c_group_open(self$uri, type = mode, ctx = soma_context(), #private$.soma_context,
                                              self$.tiledb_timestamp_range)
      }
      private$update_member_cache()
      private$update_metadata_cache()

      self
    },

    #' @description Close the SOMA object.
    #' @return The object, invisibly
    close = function() {
      if (self$is_open()) {

        for (member in private$.member_cache) {
          if (!is.null(member$object)) {
            if (member$object$is_open()) {
              member$object$close()
            } else {
            }
          }
        }

        spdl::debug("Closing {} '{}'", self$class(), self$uri)
        c_group_close(private$.tiledb_group)
        private$.mode <- NULL
        private$.tiledb_group <- NULL
        private$.group_open_timestamp <- NULL
      }
      invisible(self)
    },

    #' @description Add new member to the group. (lifecycle: maturing)
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
      spdl::debug("[TileDBGroup$set] '{}' uri {} relative {}", name, uri, relative)

      name <- name %||% basename(uri)

      private$check_open_for_write()

      #tiledb::tiledb_group_add_member(
      #  grp = private$.tiledb_group,
      #  uri = uri,
      #  relative = relative,
      #  name = name
      #)
      c_group_set(private$.tiledb_group, uri,
                  0, # -> use 'automatic' as opposed to 'relative' or 'absolute'
                  name,
                  if (inherits(object, "TileDBGroup")) "SOMAGroup" else "SOMAArray")

      private$add_cached_member(name, object)
   },

    #' @description Retrieve a group member by name. If the member isn't already
    #' open, it is opened in the same mode as the parent. (lifecycle: maturing)
    #' @param name The name of the member.
    #' @param mode Mode to open in
    #' @returns A `TileDBArray` or `TileDBGroup`.
    get = function(name) {
      stopifnot(is_scalar_character(name))
      private$check_open_for_read_or_write()

      private$fill_member_cache_if_null()
      member <- private$.member_cache[[name]]
      if (is.null(member)) {
        stop(sprintf("No member named '%s' found", name), call. = FALSE)
      }

      # Important use case:
      # * measurement$X <- SOMACollectionCreate()
      #   - At this point measurement$X is open for write
      # * measurement$X$something()
      #   - That needs to get the _actual object which is opened for write_.
      # So here if the object (maybe opened for read or write) was stored,
      # we return it. But if not (e.g. first access on read from storage)
      # then we invoke the appropriate constructor. Note: child classes
      # may override construct_member.
      obj <- if (is.null(member$object)) {
        spdl::debug("[TileDBGroup$get] construct member {} type {}", member$uri, member$type)
        obj <- private$construct_member(member$uri, member$type)
      } else {
        member$object
      }
      spdl::debug("[TileDBGroup$get] open check, mode {}", self$mode())
      if (!obj$is_open()) {
        switch(
          EXPR = (mode <- self$mode()),
          READ = obj$open(mode, internal_use_only = "allowed_use"),
          WRITE = obj$reopen(mode)
        )
      }

      private$add_cached_member(name, obj)
      obj
    },

    #' @description Remove member. (lifecycle: maturing)
    #' @param name Name of the member to remove.
    #' @export
    remove = function(name) {
      stopifnot(is_scalar_character(name))
      private$check_open_for_write()

      #tiledb::tiledb_group_remove_member(grp = private$.tiledb_group, uri = name)
      c_group_remove_member(private$.tiledb_group, name)

      # Drop member if cache has been initialized
      if (is.list(private$.member_cache)) {
        private$.member_cache[[name]] <- NULL
      }
    },

    #' @description Length in the number of members. (lifecycle: maturing)
    #' @return Scalar `integer`
    length = function() {
      private$check_open_for_read_or_write()
      private$fill_member_cache_if_null()
      length(private$.member_cache)
    },

    #' @description Retrieve the names of members. (lifecycle: maturing)
    #' @return A `character` vector of member names.
    names = function() {
      private$check_open_for_read_or_write()
      private$fill_member_cache_if_null()
      names(private$.member_cache) %||% character(length = 0L)
    },

    #' @description Retrieve a `list` of members. (lifecycle: maturing)
    to_list = function() {
      private$check_open_for_read_or_write()
      private$fill_member_cache_if_null()
      private$.member_cache
    },

    #' @description Retrieve a `data.frame` of members. (lifecycle: maturing)
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

    #' @description Retrieve metadata. (lifecycle: maturing)
    #' @param key The name of the metadata attribute to retrieve.
    #'   is not NULL.
    #' @return A list of metadata values.
    get_metadata = function(key = NULL) {
      private$check_open_for_read_or_write()

      private$fill_metadata_cache_if_null()

      spdl::debug("Retrieving metadata for {} '{}'", self$class(), self$uri)
      if (!is.null(key)) {
        val <- private$.metadata_cache[[key]]
        if (is.list(val)) val <- unlist(val)
        val
      } else {
        private$.metadata_cache
      }
    },

    #' @description Add list of metadata. (lifecycle: maturing)
    #' @param metadata Named list of metadata to add.
    #' @return NULL
    set_metadata = function(metadata) {
      stopifnot(
        "Metadata must be a named list" = is_named_list(metadata)
      )

      private$check_open_for_write()

      spdl::debug("Writing metadata to {} '{}'", self$class(), self$uri)
      dev_null <- mapply(
        FUN = c_group_put_metadata,
        key = names(metadata),
        obj = metadata,
        MoreArgs = list(xp = private$.tiledb_group),
        SIMPLIFY = FALSE
      )

      dev_null <- mapply(
        FUN = private$add_cached_metadata,
        key = names(metadata),
        val = metadata,
        SIMPLIFY = FALSE
      )
    }
  ),

  private = list(

    # @description This is a handle at the TileDB-R level
    #
    # Important implementation note:
    # * In TileDB-R there is an unopened handle obtained by tiledb::tiledb_array, which takes
    #   a URI as its argument.
    # * One may then open and close this using tiledb::tiledb_array_open (for read or write)
    #   and tiledb::tiledb_array_close, which take a tiledb_array handle as their first argument.
    #
    # However, for groups:
    # * tiledb::tiledb_group and tiledb::group_open both return an object opened for read or write.
    # * Therefore for groups we cannot imitate the behavior for arrays.
    #
    # For this reason there is a limit to how much handle-abstraction we can do in the TileDBObject
    # parent class. In particular, we cannot have a single .tiledb_object shared by both TileDBArray
    # and TileDBGroup.
    .tiledb_group = NULL,

    # This field stores the timestamp with which we opened the group, whether we used the
    # opener-supplied self$tiledb_timestamp, or defaulted to the time of opening, or neither
    # (NULL). This is the timestamp to propagate to accessed members.
    .group_open_timestamp = NULL,

    # @description List of cached group members
    # Initially NULL, once the group is created or opened, this is populated
    # with a list that's empty or contains the group members.
    .member_cache = NULL,

    # Initially NULL, once the group is created or opened, this is populated
    # with a list that's empty or contains the group metadata.
    .metadata_cache = NULL,

    ## soma_context
    .soma_context = NULL,

    # Instantiate a group member object.
    # Responsible for calling the appropriate R6 class constructor.
    construct_member = function(uri, type) {
      stopifnot(
        is_scalar_character(uri),
        is_scalar_character(type)
      )
      spdl::warn("[TileDBGroup$construct_member] uri {} type {}", uri, type)
      constructor <- switch(type,
        ARRAY     = TileDBArray$new,
        SOMAArray = TileDBArray$new,
        GROUP     = TileDBGroup$new,
        SOMAGroup = TileDBGroup$new,
        stop(sprintf("Unknown member type: %s", type), call. = FALSE)
      )
      obj <- constructor(uri, tiledbsoma_ctx = self$tiledbsoma_ctx, tiledb_timestamp = private$.group_open_timestamp,
                         platform_config = self$platform_config, internal_use_only = "allowed_use")
      obj
    },

    # ----------------------------------------------------------------
    # Important implementation note about caching:
    #
    # The caching layer is not solely a performance-enhancer. It's a necessary part of the
    # implementation.
    #
    # At the SOMA application level, the SOMA spec requires that we allow users to read array schema
    # and metadata even when the array is opened for write, and it requires that we allow users to
    # read group members and metadata even when the group is opened for write.
    #
    # At the TileDB implementation level, for arrays, we can read the array schema but not the array
    # metadata if the array is opened for write. For groups in TileDB-R, we can read the group
    # member list or the group metadata if the group is opened for read but we cannot access them
    # when the group is opened for write.

    # ----------------------------------------------------------------
    # Member caching

    # @description Retrieve all group members. (lifecycle: maturing)
    # @return A list indexed by group member names where each element is a
    # list with names: name, uri, and type.
    get_all_members_uncached_read = function(group_handle) {

      #count <- tiledb::tiledb_group_member_count(group_handle)
      count <- c_group_member_count(group_handle)
      if (count == 0) return(list())

      members <- vector(mode = "list", length = count)
      if (count == 0) return(members)

      #for (i in seq_len(count)) {
      #  members[[i]] <- setNames(
      #    # returns character vector with <type, uri, name>
      #    object = as.list(tiledb::tiledb_group_member(group_handle, i - 1L)),
      #    nm = c("type", "uri", "name")
      #  )
      #}

      # Key the list by group member name
      #names(members) <- vapply_char(members, FUN = getElement, name = "name")

      ## names list of lists
      members <- c_group_members(group_handle)

      members
    },

    fill_member_cache_if_null = function() {
      if (is.null(private$.member_cache)) {
        private$update_member_cache()
      }
    },

    update_member_cache = function() {
      spdl::debug("[TileDBGroup$updating_member_cache] class {} uri '{}'", self$class(), self$uri)

      # See notes above -- at the TileDB implementation level, we cannot read anything about the
      # group while the group is open for read, but at the SOMA application level we must support
      # this. Therefore if the group is opened for write and there is no cache populated then
      # we must open a temporary handle for read, to fill the cache.
      group_handle <- private$.tiledb_group
      if (private$.mode == "WRITE") {
        # FIXME: when we add write timestamps we should open this temp handle with tiledb_timestamp
        # too. The stopifnot is currently "unreachable" since open() stops if called with WRITE
        # mode and non-null tiledb_timestamp.
        #stopifnot("FIXME" = is.null(private$.group_open_timestamp))
        # group_handle <- tiledb::tiledb_group(self$uri, type = "READ", ctx = private$.tiledb_ctx)
        if (is.null(private$.soma_context)) private$.soma_context <- soma_context()
        spdl::debug("[TileDBGroup$updating_member_cache] re-opening {} uri '{}' ctx null {} time null {}",
                    self$class(), self$uri, is.null(private$.soma_context),
                    is.null(private$.tiledb_timestamp_range))
        group_handle <- c_group_open(self$uri, type = "READ",
                                     ctx = private$.soma_context,
                                     private$.tiledb_timestamp_range)
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
        #tiledb::tiledb_group_close(group_handle)
        c_group_close(group_handle)
      }
    },

    add_cached_member = function(name, object) {
      # We explicitly add the new member to member_cache in order to preserve the
      # original URI. Otherwise TileDB Cloud creation URIs are retrieved from
      # using tiledb_group_member() in the form tiledb://namespace/uuid. In this
      # form it's not possible to append new children, which is necessary during
      # ingestion.
      if (is.null(private$.member_cache)) {
        private$.member_cache <- list()
      }

      private$.member_cache[[name]] <- list(
        # TODO: do we really need the type here?
        # Calling tiledb::tiledb_object_type on remote storage has a cost;
        # perhaps unnecessary to incur.
        type = tiledb::tiledb_object_type(object$uri),
        uri = object$uri,
        name = name,
        object = object
      )

      # We still need to update member_cache to pick up existing members.
      # Otherwise if you open a group with existing members and add a new
      # member, the initially empty member_cache will only contain the new
      # member.
      private$update_member_cache()
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
    },

    # ----------------------------------------------------------------
    # Metadata-caching

    fill_metadata_cache_if_null = function() {
      if (is.null(private$.metadata_cache)) {
        private$update_metadata_cache()
      }
    },

    update_metadata_cache = function() {
      spdl::debug("Updating metadata cache for {} '{}'", self$class(), self$uri)

      # See notes above -- at the TileDB implementation level, we cannot read anything about the
      # group while the group is open for read, but at the SOMA application level we must support
      # this. Therefore if the group is opened for write and there is no cache populated then
      # we must open a temporary handle for read, to fill the cache.
      group_handle <- private$.tiledb_group
      if (private$.mode == "WRITE") {
        # group_handle <- tiledb::tiledb_group(self$uri, type = "READ", ctx = private$.tiledb_ctx)
        group_handle <- c_group_open(self$uri, type ="READ", ctx = private$.soma_context,
                                     private$.tiledb_timestamp_range)

      }

      #private$.metadata_cache <- tiledb::tiledb_group_get_all_metadata(group_handle)
      private$.metadata_cache <- c_group_get_metadata(group_handle)
      if (private$.mode == "WRITE") {
          #tiledb::tiledb_group_close(group_handle)
          c_group_get_metadata(group_handle)
      }

      invisible(NULL)
    },

    add_cached_metadata = function(key, value) {
      if (is.null(private$.metadata_cache)) {
        private$.metadata_cache <- list()
      }
      private$.metadata_cache[[key]] <- value
    }

  )
)
