#' The SOMA Object Base Class
#'
#' @description Base class to implement shared functionality across the
#' \code{\link{SOMAArrayBase}} and \code{\link{SOMACollectionBase}} classes.
#' (lifecycle: maturing)
#'
#' @keywords internal
#'
#' @export
#'
#' @seealso Derived classes: \code{\link{SOMAArrayBase}},
#' \code{\link{SOMACollectionBase}}
#'
SOMAObject <- R6::R6Class(
  classname = "SOMAObject",
  public = list(
    #' @description Create a new SOMA object. (lifecycle: maturing)
    #'
    #' @param uri URI for the SOMA object
    #' @param ... Ignored
    #' @param platform_config Optional platform configuration
    #' @param tiledbsoma_ctx Optional TileDB SOMA context
    #' @param tiledb_timestamp Optional timestamp (\code{\link[base]{POSIXct}})
    #' to open the object at
    #' @param soma_context A SOMA context as created by
    #' \code{\link{soma_context}()}
    #'
    initialize = function(
      uri,
      ...,
      platform_config = NULL,
      tiledbsoma_ctx = NULL,
      tiledb_timestamp = NULL,
      soma_context = NULL
    ) {
      envs <- unique(vapply(
        X = unique(sys.parents()),
        FUN = function(n) environmentName(environment(sys.function(n))),
        FUN.VALUE = character(1L)
      ))
      if (!"tiledbsoma" %in% envs) {
        msg <- private$.internal_use_only(
          "new",
          ifelse(inherits(self, "SOMACollectionBase"), "collection", "array")
        )
        stop(paste(strwrap(msg), collapse = '\n'), call. = FALSE)
      }

      # Set the URI
      if (!is.character(uri) || length(uri) != 1L || !nzchar(uri)) {
        stop("'uri' must be a single, non-empty string", call. = FALSE)
      }
      private$.uri <- uri

      # Set the platform config
      platform_config <- platform_config %||% PlatformConfig$new()
      if (!inherits(platform_config, "PlatformConfig")) {
        stop("'platform_config' must be a PlatformConfig object", call. = FALSE)
      }
      private$.platform_config <- platform_config

      # Set the context
      if (!is.null(x = tiledbsoma_ctx)) {
        if (!inherits(x = tiledbsoma_ctx, what = 'SOMATileDBContext')) {
          stop(
            "'tiledbsoma_ctx' must be a SOMATileDBContext object",
            call. = FALSE
          )
        }
        # TODO: Deprecate tiledbsoma_ctx in favor of soma_context
        # warning("'tiledbsoma_ctx' is deprecated, use 'soma_context' instead")
        # Set the old context
        private$.tiledbsoma_ctx <- tiledbsoma_ctx
        private$.tiledb_ctx <- self$tiledbsoma_ctx$context()
        # Also plumb through to the new context
        if (!is.null(soma_context)) {
          warning(
            "Both 'soma_context' and 'tiledbsoma_ctx' were provided,",
            "using 'soma_context' only"
          )
        } else {
          # why we named the parameter and function the same thing is beyond me
          soma_context <- tiledbsoma::soma_context(
            config = unlist(tiledbsoma_ctx$to_list())
          )
        }
      } else {
        tiledbsoma_ctx <- SOMATileDBContext$new()
        private$.tiledbsoma_ctx <- tiledbsoma_ctx
        private$.tiledb_ctx <- self$tiledbsoma_ctx$context()
      }

      # TODO: re-enable once new UX is worked out
      # soma_context <- soma_context %||% soma_context()
      # stopifnot(
      #   "'soma_context' must be a pointer" = inherits(x = soma_context, what = 'externalptr')
      # )
      if (is.null(soma_context)) {
        private$.soma_context <- soma_context() # FIXME via factory and paramater_config
      } else {
        private$.soma_context <- soma_context
      }

      # Set the timestamp
      if (!is.null(tiledb_timestamp)) {
        stopifnot(
          "'tiledb_timestamp' must be a single POSIXct datetime object" = inherits(tiledb_timestamp, "POSIXct") &&
            length(tiledb_timestamp) == 1L &&
            !is.na(tiledb_timestamp)
        )
        private$.tiledb_timestamp <- tiledb_timestamp
      }

      spdl::debug(
        "[SOMAObject] initialize {} with '{}' at ({})",
        self$class(),
        self$uri,
        self$tiledb_timestamp %||% "now"
      )
    },

    # NOTE: The create/open/close are necessarily specific to arrays and
    # collections; this is a bit of re-use at the ABC level

    #' @description Determine if the object is open for reading or writing
    #'
    #' @return \code{TRUE} if the object is open, otherwise \code{FALSE}
    #'
    is_open = \() self$mode() != "CLOSED",

    #' @description Print the name of the R6 class
    #'
    #' @return The name of the R6 class
    #'
    class = \() class(self)[1L],

    #' @description Get the mode of the object
    #'
    #' @return The mode of the object, one of:
    #' \itemize{
    #'  \item \dQuote{\code{CLOSED}}
    #'  \item \dQuote{\code{READ}}
    #'  \item \dQuote{\code{WRITE}}
    #'  \item \dQuote{\code{DELETE}}
    #' }
    #'
    mode = \() private$.mode %||% "CLOSED",

    #' @description Close and reopen the TileDB object in a new mode
    #'
    #' @param mode New mode to open the object in; choose from:
    #' \itemize{
    #'  \item \dQuote{\code{READ}}
    #'  \item \dQuote{\code{WRITE}}
    #'  \item \dQuote{\code{DELETE}}
    #' }
    #' @param tiledb_timestamp Optional Datetime (POSIXct) with TileDB timestamp
    #'
    #' @return Invisibly returns \code{self} opened in \code{mode}
    #'
    reopen = function(mode, tiledb_timestamp = NULL) {
      mode <- match.arg(mode, choices = c("READ", "WRITE", "DELETE"))
      stopifnot(
        "'tiledb_timestamp' must be a POSIXct datetime object" = is.null(tiledb_timestamp) ||
          (inherits(tiledb_timestamp, what = "POSIXct") && length(tiledb_timestamp) == 1L && !is.na(tiledb_timestamp))
      )
      self$close()
      private$.tiledb_timestamp <- tiledb_timestamp
      self$open(mode)
      return(invisible(self))
    },

    #' @description Check if the object exists. (lifecycle: maturing)
    #'
    #' @return \code{TRUE} if the object exists, otherwise \code{FALSE}
    #'
    exists = function() {
      expected_type <- if (self$class() == "SOMAObject") {
        c("ARRAY", "GROUP")
      } else if (inherits(self, "SOMAArrayBase")) {
        "ARRAY"
      } else if (inherits(self, "SOMACollectionBase")) {
        "GROUP"
      } else {
        stop("Unknown object type", call. = FALSE)
      }
      return(get_tiledb_object_type(self$uri, ctxxp = private$.soma_context) %in% expected_type)
    },

    #' @description Retrieve metadata. (lifecycle: maturing)
    #'
    #' @param key The name of the metadata attribute to retrieve
    #' is not NULL.
    #'
    #' @return A list of metadata values.
    #'
    get_metadata = function(key = NULL) {
      if (!(is.null(key) || (is_scalar_character(key) && nzchar(key)))) {
        stop("'key' must be a single, non-empty string", call. = FALSE)
      }
      private$.check_open()
      private$.update_metadata_cache()

      spdl::debug("Retrieving metadata for {} '{}'", self$class(), self$uri)

      if (is.null(key)) {
        return(private$.metadata_cache)
      }
      val <- private$.metadata_cache[[key]]
      if (is.list(val)) {
        val <- unlist(val)
      }
      return(val)
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

      for (i in seq_along(metadata)) {
        key <- names(metadata)[i]
        value <- metadata[[i]]
        spdl::debug(
          "[SOMAObject$set_metadata] setting key {} to {} ({})",
          key,
          value,
          class(value)
        )
        set_metadata(
          uri = self$uri,
          key = key,
          valuesxp = value,
          type = class(value),
          is_array = inherits(self, "SOMAArrayBase"),
          ctxxp = private$.soma_context,
          tsvec = self$.tiledb_timestamp_range
        )
        private$.metadata_cache[[key]] <- value
      }

      return(invisible(self))
    },

    #' @description Print-friendly representation of the object
    #'
    #' @return Invisibly returns \code{self}
    #'
    print = function() {
      cat(
        sprintf(fmt = "<%s>", self$class()),
        sprintf(fmt = "  uri: %s", self$uri),
        sep = "\n"
      )
      return(invisible(self))
    }

  ),
  active = list(
    #' @field platform_config Platform configuration
    #'
    platform_config = function(value) {
      if (!missing(x = value)) {
        private$.read_only_error("platform_config")
      }
      return(private$.platform_config)
    },

    #' @field tiledbsoma_ctx SOMATileDBContext
    #'
    tiledbsoma_ctx = function(value) {
      if (!missing(x = value)) {
        private$.read_only_error("tiledbsoma_ctx")
      }
      return(private$.tiledbsoma_ctx)
    },

    #' @field tiledb_timestamp Time that object was opened at
    #'
    tiledb_timestamp = function(value) {
      if (!missing(value)) {
        private$.read_only_error("tiledb_timestamp")
      }
      return(private$.tiledb_timestamp)
    },

    #' @field uri The URI of the TileDB object
    #'
    uri = function(value) {
      if (!missing(value)) {
        private$.read_only_error("uri")
      }
      return(private$.uri)
    },

    #' @field soma_type The SOMA object type
    #'
    soma_type = function(value) {
      if (!missing(value)) {
        private$.read_only_error("soma_type")
      }
      if (!length(private$.soma_type) || !nzchar(private$.soma_type)) {
        private$.soma_type <- self$get_metadata(SOMA_OBJECT_TYPE_METADATA_KEY)
      }
      return(private$.soma_type)
    },

    #' @field .tiledb_timestamp_range Time range for libtiledbsoma
    #'
    .tiledb_timestamp_range = function(value) {
      if (!missing(value)) {
        private$.read_only_error("tiledb_timestamp_range")
      }
      if (is.null(self$tiledb_timestamp)) {
        return(NULL)
      }
      return(c(
        as.POSIXct(0, tz = "UTC", origin = "1970-01-01"),
        self$tiledb_timestamp
      ))
    }
  ),
  private = list(

    # @field .platform_config ...
    #
    .platform_config = NULL,

    # @field .tiledbsoma_ctx ...
    #
    .tiledbsoma_ctx = NULL,

    # @field .tiledb_ctx ...
    #
    .tiledb_ctx = NULL,

    # @field .tiledb_timestamp ...
    #
    .tiledb_timestamp = NULL,

    # @field .soma_context ...
    #
    .soma_context = NULL,

    # @field .soma_type ...
    #
    .soma_type = character(1L),

    # @field .mode ...
    #
    .mode = character(1L),

    # @field .uri ...
    #
    .uri = character(1L),

    # @field .metadata_cache ...
    #
    .metadata_cache = NULL,

    # @description Create a message saying that a method is for
    # internal use only
    #
    .internal_use_only = \(method, type = c("array", "collection")) sprintf(
      fmt = "Use of the '%s()' method is for internal use only; consider using a factory method such as '%s()' instead",
      method,
      sprintf(
        fmt = "%s%s",
        switch(type, collection = "SOMACollection", "SOMADataFrame"),
        switch(method, create = "Create", "Open")
      )
    ),

    # @description Throw an error saying a field is read-only
    #
    .read_only_error = \(field) stop(
      "Field ",
      sQuote(field),
      " is read-only",
      call. = FALSE
    ),

    # @description Check that the object is open for reading
    #
    .check_open_for_read = function() {
      if (!switch(self$mode() %||% "", READ = TRUE, FALSE)) {
        stop("Item must be open for read: ", self$uri, call. = FALSE)
      }
      return(invisible(NULL))
    },

    # @description Check that the object is open for writing
    #
    .check_open_for_write = function() {
      if (!switch(self$mode() %||% "", WRITE = TRUE, FALSE)) {
        stop("Item must be open for write: ", self$uri, call. = FALSE)
      }
      return(invisible(NULL))
    },

    # @desciption Check that the object is open for delete
    .check_open_for_delete = function () {
      if (self$mode() != "DELETE") {
        stop("Item must be open for delete: ", self$uri, call. = FALSE)
      }
      return(invisible(NULL))
    },

    # @description Check that the object is open
    #
    .check_open_for_read_or_write = function() {
      if (!switch(self$mode() %||% "", READ = , WRITE = TRUE, FALSE)) {
        stop("Item must be open for read or write: ", self$uri, call. = FALSE)
      }
      return(invisible(NULL))
    },

    # @description Check that the object is open
    #
    .check_open = function() {
      if (!self$is_open()) {
        stop("Item must be open for read, write, or, delete: ", self$uri, call. = FALSE)
      }
      return(invisible(NULL))
    },

    # @description Update the metadata cache
    #
    # @param force \code{TRUE} or \code{FALSE}
    #
    # @return Invisibly returns \code{self}
    #
    .update_metadata_cache = function(force = FALSE) {
      stopifnot(isTRUE(force) || isFALSE(force))

      if (is.null(private$.metadata_cache)) {
        private$.metadata_cache <- list()
      }

      # Skip if we already have a member cache and don't want to update
      if (length(private$.metadata_cache) && !force) {
        return(invisible(NULL))
      }

      spdl::debug(
        "[SOMAObject$update_metadata_cache] updating metadata cache for {} '{}' in {}",
        self$class(),
        self$uri,
        self$mode()
      )

      private$.metadata_cache <- get_all_metadata(
        uri = self$uri,
        is_array = inherits(self, "SOMAArrayBase"),
        ctxxp = private$.soma_context
      ) %||% list()

      # Allow method chaining
      return(invisible(self))
    }
  )
)
