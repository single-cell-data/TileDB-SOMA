#' SOMA Collection Base Class
#'
#' @description Base class for objects containing persistent collection of SOMA
#' objects, mapping string keys to any SOMA object (lifecycle: maturing).
#'
#' @details Unlike \link[base]{list}s in \R, SOMA collections are unordered. The
#' first item added to a collection is not guaranteed to be the first item
#' listed in the collection. Take care when looking through the \code{$names()}
#' of a SOMA collection, or use \code{\link[base:match]{\%in\%}} to assert
#' that a specific member is present in a SOMA collection.
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
    #' @description Create a SOMA collection (lifecycle: maturing).\cr
    #' \cr
    #' \strong{Note}: \code{$create()} is considered internal and should not be
    #' called directly; use factory functions
    #' (eg. \code{\link{SOMACollectionCreate}()}) instead.
    #'
    #' @return Returns \code{self}.
    #'
    create = function() {
      private$.check_call_is_internal(
        "create",
        paste(self$class(), "Create", sep = "")
      )
      soma_debug(sprintf(
        "[SOMACollectionBase$create] Creating new %s at '%s' at %s",
        self$class(),
        self$uri,
        self$tiledb_timestamp %||% "now"
      ))
      private$.create()
      self$open("WRITE")
      return(self)
    },

    #' @description Open the SOMA collection for read or write.\cr
    #' \cr
    #' \strong{Note}: \code{$open()} is considered internal and should not be
    #' called#' directly; use factory functions
    #' (eg. \code{\link{SOMACollectionOpen}()}) instead.
    #'
    #' @param mode Mode to open the object in.
    #'
    #' @return Returns \code{self}.
    #'
    open = function(mode = c("READ", "WRITE", "DELETE")) {
      private$.check_call_is_internal(
        "open",
        paste(self$class(), "Open", sep = "")
      )
      open_mode <- match.arg(mode)
      private$.log_open_timestamp(open_mode)
      if (is.null(super$handle)) {
        private$.open_handle(open_mode, self$tiledb_timestamp)
      } else {
        soma_object_open(super$handle, open_mode, self$tiledb_timestamp)
      }
      return(self)
    },

    #' @description Close the SOMA collection.
    #'
    #' @return Invisibly returns \code{self}.
    #'
    close = function(recursive = TRUE) {
      soma_debug(sprintf(
        "[SOMAObject$close] Closing %s '%s'",
        self$class(),
        self$uri
      ))
      if (!is.null(super$handle)) {
        soma_object_close(super$handle, recursive)
      }
      return(invisible(self))
    },

    #' @description Add a new SOMA object to the collection
    #' (lifecycle: maturing).
    #'
    #' This method adds an existing SOMA object to the collection under the
    #' specified key. Replacing an existing key is not supported; attempting to
    #' add an object with a key that already exists will raise an error.
    #'
    #' @note This method is not supported for Carrara (TileDB v3) URIs. For
    #' Carrara collections, use the `add_new_*` methods instead, which create
    #' child objects at nested URIs that are automatically registered with the
    #' parent collection.
    #'
    #' @param object SOMA object.
    #' @param name The name to use for the object; defaults to the basename of
    #' \code{object$uri}.
    #' @param relative An optional logical value indicating whether the new
    #' object's URI is relative to the collection's URI. If \code{NULL} (the
    #' default), the object's URI is assumed to be relative unless it is a
    #' \code{tiledb://} URI.
    #'
    #' @return Invisibly returns \code{self}.
    #'
    set = function(object, name = NULL, relative = NULL) {
      stopifnot(
        "Only SOMA objects may be added" = inherits(object, "SOMAObject"),
        "'name' must be a single, non-empty string" = is.null(name) ||
          (is_scalar_character(name) && nzchar(name))
      )

      private$.check_open_for_write()

      # Default name to URI basename
      name <- name %||% basename(object$uri)

      # Prevent replacing existing keys with new objects
      if (name %in% self$names()) {
        stop(
          sprintf("replacing key '%s' is unsupported", name),
          call. = FALSE
        )
      }

      # Determine whether to use relative URI
      relative <- relative %||% is_relative_uri(object$uri)
      if (!(isTRUE(relative) || isFALSE(relative))) {
        stop("'relative' must be TRUE or FALSE", call. = FALSE)
      }

      # Make the URI relative before adding it
      uri <- if (relative) {
        make_uri_relative(object$uri, self$uri)
      } else {
        object$uri
      }

      soma_debug(sprintf(
        "[SOMACollectionBase$set] '%s' uri %s relative %s",
        name,
        uri,
        relative
      ))

      soma_group_set(
        xp = super$handle,
        uri = uri,
        uri_type_int = 0, # -> use 'automatic' as opposed to 'relative' or 'absolute'
        name = name,
        soma_type = ifelse(
          test = inherits(object, "SOMACollectionBase"),
          yes = "SOMAGroup",
          no = "SOMAArray"
        ),
        object$handle,
        TRUE
      )

      return(invisible(self))
    },

    #' @description Retrieve a SOMA object by name. (lifecycle: maturing)
    #'
    #' @param name The name of the object to retrieve.
    #'
    #' @return The SOMA object stored as \code{name}.
    #'
    get = function(name) {
      if (!is.character(name) || length(name) != 1L || !nzchar(name)) {
        stop("'name' must be a single, non-empty string", call. = FALSE)
      }
      private$.check_open()

      if (!c_group_has_member(super$handle, name)) {
        stop(sprintf("No member named '%s' found", name), call. = FALSE)
      }

      handle <- soma_group_get_member(super$handle, name)

      soma_debug(sprintf(
        "[SOMACollectionBase$get] open check, mode %s",
        self$mode()
      ))

      obj <- SOMAObjectWrap(
        self$members[[name]]$uri,
        tiledb_handle = handle,
        platform_config = self$platform_config,
        context = self$context,
        tiledb_timestamp = self$tiledb_timestamp,
        tiledbsoma_ctx = private$.tiledbsoma_ctx
      )

      soma_debug(sprintf(
        "[SOMACollectionBase$get] open check, mode %s",
        self$mode()
      ))
      if (!obj$is_open()) {
        switch(
          EXPR = (mode <- self$mode()),
          READ = obj$open(mode),
          WRITE = obj$reopen(mode)
        )
      }

      return(obj)
    },

    #' @description Remove member. (lifecycle: maturing)
    #'
    #' @param name Name of the member to remove.
    #'
    #' @return Invisibly returns \code{self}
    #'
    remove = function(name) {
      if (self$mode() == "WRITE") {
        .Deprecated(
          msg = sprintf(
            "Removing a member in %s mode is deprecated. Collection should be opened in %s mode.",
            sQuote("WRITE"),
            sQuote("DELETE")
          )
        )
      } else if (self$mode() != "DELETE") {
        stop(
          "SOMA object is not opened in 'delete' mode; cannot remove member.",
          call. = FALSE
        )
      }

      if (!is.character(name) || length(name) != 1L || !nzchar(name)) {
        stop("'name' must be a single, non-empty string", call. = FALSE)
      }
      private$.check_open()

      soma_group_remove_member(super$handle, name)

      return(invisible(self))
    },

    #' @description Get the number of members in this collection
    #' (lifecycle: maturing).
    #'
    #' @return The number of members in this collection.
    #'
    length = function() {
      private$.check_open()
      return(length(self$members))
    },

    #' @description Retrieve the names of members (lifecycle: maturing).
    #'
    #' Note: the order of the member names is not guaranteed.
    #'
    #' @return A character vector of member names.
    #'
    names = function() {
      private$.check_open()
      return(names(self$members) %||% character(length = 0L))
    },

    #' @description Add a new SOMA collection to this collection
    #'  (lifecycle: maturing).
    #'
    #' @param object SOMA collection object.
    #' @param key The key to be added.
    #'
    #' @section Carrara (TileDB v3) behavior:
    #'
    #' For Carrara URIs, child objects created at nested URIs are automatically
    #' added to the parent collection. Calling this method on an already-
    #' registered child is a **no-op** for backward compatibility.
    #'
    #' @return Returns \code{object}.
    #'
    add_new_collection = function(object, key) {
      if (!inherits(object, "SOMACollectionBase")) {
        stop("'object' must be a SOMA collection", call. = FALSE)
      }

      private$.set_element(object, key)
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
    #' @return Returns the newly created data frame stored at \code{key}.
    #'
    add_new_dataframe = function(
      key,
      schema,
      index_column_names,
      domain,
      platform_config = NULL
    ) {
      if (key %in% self$names()) {
        stop(sprintf("Member '%s' already exists", key), call. = FALSE)
      }

      sdf <- SOMADataFrameCreate(
        uri = file_path(self$uri, key),
        schema = schema,
        index_column_names = index_column_names,
        domain = domain,
        platform_config = platform_config %||% self$platform_config,
        context = self$context,
        tiledb_timestamp = self$tiledb_timestamp # Cached value from $new()/SOMACollectionOpen
      )
      private$.set_element(sdf, key)
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
    #' @return Returns the newly-created array stored at \code{key}.
    #'
    add_new_dense_ndarray = function(key, type, shape, platform_config = NULL) {
      if (key %in% self$names()) {
        stop(sprintf("Member '%s' already exists", key), call. = FALSE)
      }

      ndarr <- SOMADenseNDArrayCreate(
        uri = file_path(self$uri, key),
        type = type,
        shape = shape,
        platform_config = platform_config %||% self$platform_config,
        context = self$context,
        tiledb_timestamp = self$tiledb_timestamp
      )
      private$.set_element(ndarr, key)
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
    #' @return Returns the newly-created array stored at \code{key}.
    #'
    add_new_sparse_ndarray = function(
      key,
      type,
      shape,
      platform_config = NULL
    ) {
      if (key %in% self$names()) {
        stop(sprintf("Member '%s' already exists", key), call. = FALSE)
      }

      ndarr <- SOMASparseNDArrayCreate(
        uri = file_path(self$uri, key),
        type = type,
        shape = shape,
        platform_config = platform_config %||% self$platform_config,
        context = self$context,
        tiledb_timestamp = self$tiledb_timestamp # Cached value from $new()/SOMACollectionOpen
      )
      private$.set_element(ndarr, key)
      return(ndarr)
    },

    #' @description Print-friendly representation of the object.
    #'
    #' @return Invisibly returns \code{self}.
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
    #' @field members A list with the members of this collection.
    #'
    members = function(value) {
      if (!missing(value)) {
        private$.read_only_error("members")
      }
      return(soma_group_get_members(super$handle))
    }
  ),
  private = list(
    # @field .group_open_timestamp ...
    #
    .group_open_timestamp = NULL,

    # @description Internal method to add a newly-created element to a
    # collection. For Carrara (v3) URIs, children are auto-registered when
    # created at a nested URI, so we only update the cache. For v2 URIs, we
    # also register the member with the TileDB group.
    #
    # @param object A SOMA object
    # @param name The key for the object
    #
    # @return Invisibly returns self
    #
    .set_element = function(object, name) {
      if (self$context$is_tiledbv3(self$uri)) {
        # Carrara requires member name to match URI basename
        if (basename(object$uri) != name) {
          stop(
            sprintf(
              paste(
                "Member name `%s` must match the final segment of the URI",
                "(`%s`) for Carrara collections."
              ),
              name,
              basename(object$uri)
            ),
            call. = FALSE
          )
        }
        # Carrara: children are auto-registered, just update cache
      }

      # v2: register with TileDB group
      self$set(object, name)
      
      return(invisible(self))
    },

    # Internal method called by SOMA Measurement/Experiment's active bindings
    # to retrieve or set one of the pre-defined SOMA fields (e.g., obs, X, etc). (lifecycle: maturing)
    # @param value the optional argument passed to the active binding.
    # @param name the name of the field to retrieve or set.
    # @param expected_class the expected class of the value to set.
    get_or_set_soma_field = function(value, name, expected_class) {
      private$.check_open()

      if (missing(value)) {
        return(self$get(name))
      }

      stopifnot(
        "Must define 'name' of the field to set" = !missing(name),
        "Must define the field's 'expected_class'" = !missing(expected_class)
      )

      if (!inherits(value, expected_class)) {
        stop(sprintf("%s must be a '%s'", name, expected_class), call. = FALSE)
      }
      self$set(value, name = name)
      return(invisible(self))
    }
  )
)
