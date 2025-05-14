#' SOMA Collection Base Class
#'
#' @description Base class for objects containing persistent collection of SOMA
#' objects, mapping string keys to any SOMA object (lifecycle: maturing).
#'
#' @keywords internal
#'
#' @export
#'
#' @seealso Derived classes: \code{\link{SOMACollection}},
#' \code{\link{SOMAMeasurement}},
#' \code{\link{SOMAExperiment}}
#'
SOMACollectionBase <- R6::R6Class(
  classname = "SOMACollectionBase",
  inherit = SOMAObject,
  public = list(

    #' @description Create a new \code{\link{SOMACollection}}
    #' (lifecycle: maturing).
    #'
    #' @param uri URI of the TileDB group.
    #' @param platform_config Optional storage-engine specific configuration.
    #' @param tiledbsoma_ctx optional SOMATileDBContext.
    #' @param tiledb_timestamp Optional Datetime (POSIXct) for TileDB timestamp.
    #' @param internal_use_only Character value to signal this is a 'permitted'
    #' call, as \code{new()} is considered internal and should not be called
    #' directly.
    #'
    initialize = function(
      uri,
      platform_config = NULL,
      tiledbsoma_ctx = NULL,
      tiledb_timestamp = NULL,
      internal_use_only = NULL
    ) {
      return(super$initialize(
        uri = uri,
        platform_config = platform_config,
        tiledbsoma_ctx = tiledbsoma_ctx,
        tiledb_timestamp = tiledb_timestamp,
        internal_use_only = internal_use_only
      ))
    },

    #' @description Add a new SOMA object to the collection
    #' (lifecycle: maturing).
    #'
    #' @param internal_use_only Character value to signal this is a 'permitted'
    #' call, as \code{create()} is considered internal and should not be called
    #' directly.
    #'
    #' @return Returns \code{self}
    #'
    create = function(internal_use_only = NULL) {
      if (is.null(internal_use_only) || internal_use_only != "allowed_use") {
        stop(paste(
          "Use of the create() method is for internal use only. Consider using a",
          "factory method as e.g. 'SOMACollectionCreate()'."
        ), call. = FALSE)
      }
      super$create(internal_use_only = internal_use_only)

      # Root SOMA objects include a `dataset_type` entry to allow the
      # TileDB Cloud UI to detect that they are SOMA datasets.
      metadata <- switch(
        self$class(),
        SOMAExperiment = list(dataset_type = "soma"),
        list()
      )

      # Key constraint on the TileDB Group API:
      # * tiledb_group_create is synchronous: new items appear on disk once that
      #   function has returned.
      # * When we have an open handle for write and we do tiledb_group_put_metadata,
      #   metadata items do _not_ appear on disk until tiledb_group_close has been
      #   called.
      # To avoid confusion, we pay the price of a reopen here.
      #
      # TODO: this feels like it could be resolved through improved caching.
      self$open("WRITE")
      private$write_object_type_metadata(metadata)
      self$reopen(mode = "WRITE")

      return(self)
    },

    #' @description Open the SOMA collection for read or write
    #'
    #' @param mode Mode to open the object in
    #' @param ... Ignored
    #'
    #' @return \code{self}
    #'
    #' @note \code{open()} is considered internal and should not be called
    #' directly; use factory functions (eg. \code{\link{SOMACollectionOpen}()})
    #' instead
    #'
    open = function(mode = c("READ", "WRITE"), ...) {
      envs <- unique(vapply(
        X = unique(sys.parents()),
        FUN = function(n) environmentName(environment(sys.function(n))),
        FUN.VALUE = character(1L)
      ))
      if (sys.parent()) {
        if (inherits(environment(sys.function(sys.parent()))$self, what = "SOMAObject")) {
          envs <- union(envs, "tiledbsoma")
        }
      }
      if (!"tiledbsoma" %in% envs) {
        stop(
          paste(strwrap(private$.internal_use_only("open", "collection")), collapse = '\n'),
          call. = FALSE
        )
      }

      # Set the mode of the collection
      private$.mode <- match.arg(mode)

      # Set the group timestamp
      # In READ mode, if the opener supplied no timestamp then we default to the time of
      # opening, providing a temporal snapshot of all group members.
      private$.group_open_timestamp <- if (self$mode() == "READ" && is.null(self$tiledbtimestamp)) {
        Sys.time()
      } else {
        self$tiledb_timestamp
      }

      private$.tiledb_group <- c_group_open(
        uri = self$uri,
        type = self$mode(),
        ctxxp = private$.soma_context,
        timestamp = self$.tiledb_timestamp_range
      )

      private$.update_member_cache(TRUE)

      return(self)
    },

    #' @description Close the SOMA collection
    #'
    #' @return Invisibly returns \code{self}
    #'
    close = function() {
      if (self$is_open()) {
        for (member in private$.member_cache) {
          if (!is.null(member$object) && member$object$is_open()) {
            member$object$close()
          }
        }

        spdl::debug("Closing {} '{}'", self$class(), self$uri)
        c_group_close(private$.tiledb_group)
        private$.mode <- NULL
        private$.tiledb_group <- NULL
        private$.group_open_timestamp <- NULL
      }

      return(invisible(self))
    },

    #' @description Add a new SOMA object to the collection
    #' (lifecycle: maturing).
    #'
    #' @param object SOMA object.
    #' @param name The name to use for the object; defaults to the basename of
    #' \code{object$uri}
    #' @param relative An optional logical value indicating whether the new
    #' object's URI is relative to the collection's URI. If \code{NULL} (the
    #' default), the object's URI is assumed to be relative unless it is a
    #' \code{tiledb://} URI.
    #'
    set = function(object, name = NULL, relative = NULL) {
      stopifnot(
        "Only SOMA objects may be added" = inherits(object, "SOMAObject"),
        "'name' must be a single, non-empty string" = is.null(name) ||
          (is_scalar_character(name) && nzchar(name))
      )
      relative <- relative %||% !startsWith(object$uri, "tiledb://")
      if (!(isTRUE(relative) || isFALSE(relative))) {
        stop("'relative' must be TRUE or FALSE", call. = FALSE)
      }

      private$.check_open_for_write()

      # Make the URI relative before adding it
      uri <- if (relative) make_uri_relative(object$uri, self$uri) else object$uri
      name <- name %||% basename(uri)
      spdl::debug(
        "[SOMACollectionBase$set] '{}' uri {} relative {}",
        name,
        uri,
        relative
      )

      c_group_set(
        xp = private$.tiledb_group,
        uri = uri,
        uri_type_int = 0, # -> use 'automatic' as opposed to 'relative' or 'absolute'
        name = name,
        soma_type = ifelse(
          test = inherits(object, "SOMACollectionBase"),
          yes = "SOMAGroup",
          no = "SOMAArray"
        )
      )

      private$.add_cache_member(name, object)
      return(invisible(self))
    },

    #' @description Retrieve a SOMA object by name. (lifecycle: maturing)
    #'
    #' @param name The name of the object to retrieve.
    #'
    #' @return The SOMA object stored as \code{name}.
    #'
    get = function(name) {
      # browser()
      if (!is.character(name) || length(name) != 1L || !nzchar(name)) {
        stop("'name' must be a single, non-empty string", call. = FALSE)
      }
      private$.check_open_for_read_or_write()
      private$.update_member_cache()

      if (is.null(member <- self$members[[name]])) {
        stop(sprintf("No member named '%s' found", name), call. = FALSE)
      }

      obj <- if (is.null(member$object)) {
        spdl::debug(
          "[SOMACollectionBase$get] construct member {} type {}",
          member$uri,
          member$type
        )
        private$construct_member(member$uri, member$type)
      } else {
        member$object
      }

      spdl::debug("[SOMACollectionBase$get] open check, mode {}", self$mode())
      if (!obj$is_open()) {
        switch(
          EXPR = (mode <- self$mode()),
          READ = obj$open(mode),
          WRITE = obj$reopen(mode)
        )
      }

      private$.add_cache_member(name, obj)
      return(obj)
    },

    #' @description Remove member. (lifecycle: maturing)
    #'
    #' @param name Name of the member to remove.
    #'
    #' @return Invisibly returns \code{self}
    #'
    remove = function(name) {
      if (!is.character(name) || length(name) != 1L || !nzchar(name)) {
        stop("'name' must be a single, non-empty string", call. = FALSE)
      }
      private$.check_open_for_read_or_write()

      c_group_remove_member(private$.tiledb_group, name)

      # Drop member if cache has been initialized
      if (length(self$members)) {
        private$.member_cache[[name]] <- NULL
      }

      return(invisible(self))
    },

    #' @description Get the number of members in this collection (lifecycle: maturing)
    #'
    #' @return The number of members in this collection
    #'
    length = function() {
      private$.check_open_for_read_or_write()
      private$.update_member_cache()
      return(length(self$members))
    },

    #' @description Retrieve the names of members. (lifecycle: maturing)
    #'
    #' @return A character vector of member names.
    #'
    names = function() {
      private$.check_open_for_read_or_write()
      private$.update_member_cache()
      names(self$members) %||% character(length = 0L)
    },

    #' @description Add list of metadata. (lifecycle: maturing)
    #'
    #' @param metadata Named list of metadata to add
    #'
    #' @return Invisibly returns \code{self}
    #'
    set_metadata = function(metadata) {
      stopifnot(
        "Metadata must be a named list" = is_named_list(metadata)
      )

      private$.check_open_for_write()
      private$.update_metadata_cache()

      spdl::debug("Writing metadata to {} '{}'", self$class(), self$uri)
      for (i in seq_along(metadata)) {
        key <- names(metadata)[i]
        obj <- metadata[[i]]
        c_group_put_metadata(xp = private$.tiledb_group, key = key, obj = obj)
        private$.metadata_cache[[key]] <- obj
      }

      return(invisible(self))
    },

    #' @description Add a new SOMA collection to this collection
    #'  (lifecycle: maturing).
    #'
    #' @param object SOMA collection object.
    #' @param key The key to be added.
    #'
    #' @return Returns \code{object}
    #'
    add_new_collection = function(object, key) {
      if (!inherits(object, "SOMACollectionBase")) {
        stop("'object' must be a SOMA collection", call. = FALSE)
      }
      self$set(object, key)
      return(object)
    },

    #' @description Add a new SOMA data frame to this collection
    #' (lifecycle: maturing).
    #'
    #' @param key The key to be added.
    #' @param schema Arrow schema argument passed on to
    #' \code{SOMADataFrame$create()}.
    #' @param index_column_names Index column names passed on to
    #' \code{SOMADataFrame$create()}.
    #' @param domain As in \code{\link{SOMADataFrameCreate}}.
    #' @template param-platform-config
    #'
    #' @return Returns the newly created data frame stored at \code{key}
    #'
    add_new_dataframe = function(
      key,
      schema,
      index_column_names,
      domain,
      platform_config = NULL
    ) {
      sdf <- SOMADataFrameCreate(
        uri = file_path(self$uri, key),
        schema = schema,
        index_column_names = index_column_names,
        domain = domain,
        platform_config = platform_config %||% self$platform_config,
        tiledbsoma_ctx = self$tiledbsoma_ctx,
        tiledb_timestamp = self$tiledb_timestamp # Cached value from $new()/SOMACollectionOpen
      )
      self$set(sdf, key)
      return(sdf)
    },

    #' @description Add a new SOMA DenseNdArray to this collection
    #' (lifecycle: maturing).
    #'
    #' @param key The key to be added.
    #' @param type An \link[arrow:data-type]{Arrow type} defining the
    #' type of each element in the array.
    #' @param shape a vector of integers defining the shape of the array.
    #' @template param-platform-config
    #'
    #' @return Returns the newly-created array stored at \code{key}
    #'
    add_new_dense_ndarray = function(key, type, shape, platform_config = NULL) {
      ndarr <- SOMADenseNDArrayCreate(
        uri = file_path(self$uri, key),
        type = type,
        shape = shape,
        platform_config = platform_config %||% self$platform_config,
        tiledbsoma_ctx = self$tiledbsoma_ctx,
        tiledb_timestamp = self$tiledb_timestamp
      )
      self$set(ndarr, key)
      return(ndarr)
    },

    #' @description Add a new SOMA SparseNdArray to this collection
    #' (lifecycle: maturing).
    #'
    #' @param key The key to be added.
    #' @param type An \link[arrow:data-type]{Arrow type} defining the
    #' type of each element in the array.
    #' @param shape a vector of integers defining the shape of the array.
    #' @template param-platform-config
    #'
    #' @return Returns the newly-created array stored at \code{key}
    #'
    add_new_sparse_ndarray = function(key, type, shape, platform_config = NULL) {
      ndarr <- SOMASparseNDArrayCreate(
        uri = file_path(self$uri, key),
        type = type,
        shape = shape,
        platform_config = platform_config %||% self$platform_config,
        tiledbsoma_ctx = self$tiledbsoma_ctx,
        tiledb_timestamp = self$tiledb_timestamp # Cached value from $new()/SOMACollectionOpen
      )
      self$set(ndarr, key)
      return(ndarr)
    },

    #' @description Print-friendly representation of the object
    #'
    #' @return Invisibly returns \code{self}
    #'
    print = function() {
      super$print()
      if (self$exists()) {
        if (self$is_open()) {
          if (self$length()) {
            df <- vector(mode = "list", length = 3L)
            names(df) <- c("name", "uri", "type")
            for (i in names(df)) {
              df[[i]] <- vapply(
                X = self$members,
                FUN = `[[`,
                FUN.VALUE = character(length = 1L),
                i,
                USE.NAMES = FALSE
              )
            }
            members <- paste0(
              df$name,
              ifelse(is_remote_uri(df$uri), yes = "*", no = "")
            )
            members <- split(members, df$type)
            if (!is.null(members$ARRAY)) {
              cat("  arrays:", string_collapse(sort(members$ARRAY)), "\n")
            }
            if (!is.null(members$GROUP)) {
              cat("  groups:", string_collapse(sort(members$GROUP)), "\n")
            }
          }
        } else {
          cat("  closed\n")
        }
      }
      return(invisible(self))
    }
  ),
  active = list(
    #' @field members A list with the members of this collection
    #'
    members = function(value) {
      if (!missing(value)) {
        private$.read_only_error("members")
      }
      private$.update_member_cache()
      return(private$.member_cache %||% NULL)
    }
  ),
  private = list(

    # @field .tiledb_group ...
    #
    .tiledb_group = NULL,

    # @field .group_open_timestamp ...
    #
    .group_open_timestamp = NULL,

    # @field .member_cache ...
    #
    .member_cache = NULL,

    # @description Update the member cache
    #
    # @param force \code{TRUE} or \code{FALSE}
    #
    # @return Invisibly returns \code{self}
    #
    .update_member_cache = function(force = FALSE) {
      stopifnot(isTRUE(force) || isFALSE(force))

      # Skip if we already have a member cache and don't want to update
      if (length(private$.member_cache) && !force) {
        return(invisible(NULL))
      }

      spdl::debug(
        "[SOMACollectionBase$updating_member_cache] class {} uri '{}'",
        self$class(),
        self$uri
      )

      # Get a read-handle for the group
      handle <- if (self$mode() == "WRITE") {
        spdl::debug(
          "[SOMACollectionBase$updating_member_cache] re-opening {} uri '{}' ctx null {} time null {}",
          self$class(),
          self$uri,
          is.null(private$.soma_context),
          is.null(self$.tiledb_timestamp_range)
        )
        c_group_open(
          uri = self$uri,
          type = "READ",
          ctxxp = private$.soma_context %||% soma_context(),
          timestamp = self$.tiledb_timestamp_range
        )
      } else {
        private$.tiledb_group
      }

      # Get the group members
      members <- if (!c_group_member_count(xp = handle)) {
        vector(mode = "list", length = 0L)
      } else {
        c_group_members(xp = handle)
      }

      # Don't clobber existing cache
      if (!length(private$.member_cache)) {
        private$.member_cache <- members
      } else {
        members <- members[setdiff(names(members), names(private$.member_cache))]
        private$.member_cache <- utils::modifyList(private$.member_cache, members)
      }

      # Close the read-handle if the group is open for writing
      if (self$mode() == "WRITE") {
        c_group_close(xp = handle)
      }

      # Allow method chaining
      return(invisible(self))
    },

    # @description Explicitly add a member to the cache in order to preserve
    # original URIs; otherwise, TileDB-Cloud creation URIs are retrieved which
    # prevent appending children
    #
    # @param name Name of object
    # @param object A SOMA object
    #
    # @return Invisibly returns \code{self}
    #
    .add_cache_member = function(name, object) {
      stopifnot(
        is.character(name) && length(name) == 1L && nzchar(name),
        inherits(object, "SOMAObject")
      )
      if (is.null(private$.member_cache)) {
        private$.member_cache <- list()
      }

      private$.member_cache[[name]] <- list(
        # TODO: do we really need the type here?
        # Calling `get_tiledb_object_type()` on remote storage has a cost;
        # perhaps unnecessary to incur.
        type = get_tiledb_object_type(object$uri, private$.soma_context),
        uri = object$uri,
        name = name,
        object = object
      )
      private$.update_member_cache(force = TRUE)
      return(invisible(self))
    },

    # @description Update the metadata cache
    #
    # @param force \code{TRUE} or \code{FALSE}
    #
    # @return Invisibly returns \code{self}
    #
    .update_metadata_cache = function(force = FALSE) {
      stopifnot(isTRUE(force) || isFALSE(force))

      # Skip if we already have a member cache and don't want to update
      if (length(private$.metadata_cache) && !force) {
        return(invisible(NULL))
      }

      spdl::debug("Updating metadata cache for {} '{}'", self$class(), self$uri)

      # Get a read-handle for the group
      handle <- if (self$mode() == "WRITE") {
        spdl::debug(
          "[SOMACollectionBase$updating_member_cache] re-opening {} uri '{}' ctx null {} time null {}",
          self$class(),
          self$uri,
          is.null(private$.soma_context),
          is.null(self$.tiledb_timestamp_range)
        )
        c_group_open(
          uri = self$uri,
          type = "READ",
          ctxxp = private$.soma_context %||% soma_context(),
          timestamp = self$.tiledb_timestamp_range
        )
      } else {
        private$.tiledb_group
      }

      # Get the group metadata
      private$.metadata_cache <- c_group_get_metadata(xp = handle) %||% list()

      # Close the read-handle if the group is open for writing
      if (self$mode() == "WRITE") {
        c_group_close(xp = handle)
      }

      # Allow method chaining
      return(invisible(self))
    },

    # Add standard SOMA metadata values to the object with the option to include
    # additional metadata.
    write_object_type_metadata = function(metadata = list()) {
      stopifnot(is.list(metadata))
      ## these entriess are now written by libtiledbsoma
      ##   metadata[[SOMA_OBJECT_TYPE_METADATA_KEY]] <- self$class()
      ##   metadata[[SOMA_ENCODING_VERSION_METADATA_KEY]] <- SOMA_ENCODING_VERSION
      if (length(metadata) > 0) self$set_metadata(metadata)
    },

    # Instantiate a soma member object.
    # Responsible for calling the appropriate R6 class constructor.
    construct_member = function(uri, type) {
      stopifnot(
        is_scalar_character(uri),
        is_scalar_character(type)
      )
      spdl::debug(
        "[SOMACollectionBase$construct_member] entered, uri {} type {}",
        uri,
        type
      )

      array <- switch(
        EXPR = type,
        ARRAY = ,
        SOMAArray = TRUE,
        GROUP = ,
        SOMAGroup = FALSE,
        stop("Unknown member SOMA type: ", type, call. = FALSE)
      )
      metadata <- get_all_metadata(
        uri,
        is_array = array,
        ctxxp = private$.soma_context
      )
      soma_type <- metadata$soma_object_type
      if (is.null(soma_type)) {
        stop("SOMA object type metadata is missing; cannot construct", call. = FALSE)
      }
      spdl::debug(
        "[SOMACollectionBase$construct_member] Instantiating {} object at: '{}'",
        soma_type,
        uri
      )

      fxn <- tryCatch(
        base::get(
          sprintf("%sOpen", soma_type),
          envir = getNamespace("tiledbsoma"),
          mode = "function"
        ),
        error = \() stop("Unknown member SOMA type: ", soma_type, call. = FALSE)
      )

      return(fxn(
        uri,
        mode = self$mode(),
        platform_config = self$platform_config,
        tiledbsoma_ctx = self$tiledbsoma_ctx,
        tiledb_timestamp = self$tiledb_timestamp
      ))
    },

    # Internal method called by SOMA Measurement/Experiment's active bindings
    # to retrieve or set one of the pre-defined SOMA fields (e.g., obs, X, etc). (lifecycle: maturing)
    # @param value the optional argument passed to the active binding.
    # @param name the name of the field to retrieve or set.
    # @param expected_class the expected class of the value to set.
    get_or_set_soma_field = function(value, name, expected_class) {
      private$.check_open_for_read_or_write()

      if (missing(value)) {
        return(self$get(name))
      }

      stopifnot(
        "Must define 'name' of the field to set" = !missing(name),
        "Must define the field's 'expected_class'" = !missing(expected_class)
      )

      if (!inherits(value, expected_class)) {
        stop(
          sprintf("%s must be a '%s'", name, expected_class),
          call. = FALSE
        )
      }
      self$set(value, name = name)
      return(invisible(self))
    }
  )
)
