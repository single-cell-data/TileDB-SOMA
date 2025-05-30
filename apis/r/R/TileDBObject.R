#' TileDB Object Base Class
#'
#' @description Virtual base class to implement shared functionality across the
#' TileDBArray and TileDBGroup classes (lifecycle: maturing).
#'
#' @export
#'
#' @seealso Derived classes: \code{\link{TileDBArray}},
#' \code{\link{TileDBGroup}}
#'
TileDBObject <- R6::R6Class(
  classname = "TileDBObject",
  public = list(
    #' @description Create a new TileDB object (lifecycle: maturing).
    #'
    #' @param uri URI for the TileDB object.
    #' @param platform_config Optional platform configuration.
    #' @param tiledbsoma_ctx Optional SOMATileDBContext.
    #' @param tiledb_timestamp Optional Datetime (POSIXct) with TileDB
    #' timestamp.
    #' @param internal_use_only Character value to signal this is a 'permitted'
    #' call, as \code{new()} is considered internal and should not be called
    #' directly.
    #' @param soma_context A SOMA context as created by
    #' \code{\link{soma_context}()}.
    #'
    initialize = function(
      uri,
      platform_config = NULL,
      tiledbsoma_ctx = NULL,
      tiledb_timestamp = NULL,
      internal_use_only = NULL,
      soma_context = NULL
    ) {
      if (is.null(internal_use_only) || internal_use_only != "allowed_use") {
        stop(paste(
          "Use of the new() method is for internal use only. Consider using a",
          "factory method as e.g. 'SOMADataFrameOpen()'."
        ), call. = FALSE)
      }
      if (missing(uri)) stop("Must specify a `uri`", call. = FALSE)
      # private$tiledb_uri <- TileDBURI$new(uri)
      private$tiledb_uri <- uri

      # Set platform config
      platform_config <- platform_config %||% PlatformConfig$new()
      if (!inherits(platform_config, "PlatformConfig")) {
        stop("'platform_config' must be a PlatformConfig object", call. = FALSE)
      }
      private$.tiledb_platform_config <- platform_config

      # Set context
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

      if (!is.null(tiledb_timestamp)) {
        stopifnot(
          "'tiledb_timestamp' must be a single POSIXct datetime object" = inherits(tiledb_timestamp, "POSIXct") &&
            length(tiledb_timestamp) == 1L &&
            !is.na(tiledb_timestamp)
        )
        private$.tiledb_timestamp <- tiledb_timestamp
      }

      spdl::debug(
        "[TileDBObject] initialize {} with '{}' at ({})", self$class(), self$uri,
        self$tiledb_timestamp %||% "now"
      )
    },

    #' @description Print the name of the R6 class.
    #'
    class = function() {
      class(self)[1]
    },

    # The create/open/close are necessarily specific to TileDBArray/TileDBGroup.
    # This is a bit of re-use at the TileDBObject level.
    #' @description Determine if the object is open for reading or writing.
    #'
    #' @return \code{TRUE} if the object is open, otherwise \code{FALSE}.
    #'
    is_open = function() {
      return(self$mode() != "CLOSED")
    },

    # TODO: make this an active
    #' @description Get the mode of the object.
    #'
    #' @return If the object is closed, returns \dQuote{\code{CLOSED}};
    #' otherwise returns the mode (eg. \dQuote{\code{READ}}) of the object.
    #'
    mode = function() {
      if (is.null(private$.mode)) {
        "CLOSED"
      } else {
        private$.mode
      }
    },

    #' @description Close and reopen the TileDB object in a new mode.
    #'
    #' @param mode New mode to open the object in; choose from:
    #' \itemize{
    #'  \item \dQuote{\code{READ}}
    #'  \item \dQuote{\code{WRITE}}
    #' }
    #' @param tiledb_timestamp Optional Datetime (POSIXct) with TileDB
    #' timestamp.
    #'
    #' @return Invisibly returns \code{self} opened in \code{mode}.
    #'
    reopen = function(mode, tiledb_timestamp = NULL) {
      mode <- match.arg(mode, choices = c("READ", "WRITE"))
      stopifnot(
        "'tiledb_timestamp' must be a POSIXct datetime object" = is.null(tiledb_timestamp) ||
          (inherits(tiledb_timestamp, what = "POSIXct") && length(tiledb_timestamp) == 1L && !is.na(tiledb_timestamp))
      )
      self$close()
      private$.tiledb_timestamp <- tiledb_timestamp
      self$open(mode, internal_use_only = "allowed_use")
      return(invisible(self))
    },

    #' @description Print-friendly representation of the object.
    #'
    print = function() {
      cat(glue::glue("<{self$class()}>"), sep = "\n")
      cat("  uri:", self$uri, "\n")
    },

    #' @description Check if the object exists (lifecycle: maturing).
    #'
    #' @return \code{TRUE} if the object exists, otherwise \code{FALSE}.
    #'
    exists = function() {
      if (self$class() == "TileDBObject") {
        expected_type <- c("ARRAY", "GROUP")
      } else if (inherits(self, "TileDBArray")) {
        expected_type <- "ARRAY"
      } else if (inherits(self, "TileDBGroup")) {
        expected_type <- "GROUP"
      } else {
        stop("Unknown object type", call. = FALSE)
      }
      get_tiledb_object_type(self$uri, ctx = soma_context()) %in% expected_type
    }
  ),
  active = list(
    #' @field platform_config Platform configuration.
    #'
    platform_config = function(value) {
      if (!missing(x = value)) {
        stop("'platform_config' is a read-only field", call. = FALSE)
      }
      return(private$.tiledb_platform_config)
    },

    #' @field tiledbsoma_ctx SOMATileDBContext.
    #'
    tiledbsoma_ctx = function(value) {
      if (!missing(x = value)) {
        stop("'tiledbsoma_ctx' is a read-only field", call. = FALSE)
      }
      return(private$.tiledbsoma_ctx)
    },

    #' @field tiledb_timestamp Time that object was opened at.
    #'
    tiledb_timestamp = function(value) {
      if (!missing(value)) {
        private$.read_only_error("tiledb_timestamp")
      }
      return(private$.tiledb_timestamp)
    },

    #' @field uri The URI of the TileDB object.
    #'
    uri = function(value) {
      if (missing(value)) {
        return(private$tiledb_uri)
      }
      stop(sprintf("'%s' is a read-only field.", "uri"), call. = FALSE)
    },

    #' @field .tiledb_timestamp_range Time range for libtiledbsoma.
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
    # Pro tip: in R6 we can't set these to anything other than NULL here, even
    # if we want to.  If you want them defaulted to anything other than NULL,
    # leave them NULL here and set the defaults in the constructor.

    # Set by TileDBArray and TileDBGroup, as stateful handles have incompatible
    # semanticswe can't completely abstract here in this parent class
    #
    # Semantics:
    # * "READ" when opened for read
    # * "WRITE" when opened for write
    # * NULL when never opened, or when closed.
    # * In particular, an is-open predicate can be reliably implemented by
    #   checking if .mode is non-null.
    .mode = NULL,

    ## soma_context
    .soma_context = NULL,

    # @description Contains TileDBURI object
    tiledb_uri = NULL,

    # Internal platform config
    .tiledb_platform_config = NULL,

    # Opener-supplied POSIXct timestamp, if any. TileDBArray and TileDBGroup
    # are each responsible for making this effective, since the methods
    # differ slightly.
    .tiledb_timestamp = NULL,

    # Internal context
    .tiledbsoma_ctx = NULL,
    .tiledb_ctx = NULL,
    .read_only_error = function(field) {
      stop("Field ", sQuote(field), " is read-only", call. = FALSE)
    },
    is_open_for_read = function() {
      # Pro-tip: it's not enough to check $private.mode != "READ", since logical(0) isn't
      # the same as FALSE
      if (is.null(private$.mode)) {
        FALSE
      } else if (private$.mode != "READ") {
        FALSE
      } else {
        TRUE
      }
    },
    is_open_for_write = function() {
      if (is.null(private$.mode)) {
        FALSE
      } else if (private$.mode != "WRITE") {
        FALSE
      } else {
        TRUE
      }
    },

    # Per the spec, invoking user-level read requires open for read mode.
    check_open_for_read = function() {
      if (!private$is_open_for_read()) {
        stop(paste("Item must be open for read:", self$uri))
      }
    },

    # Per the spec, invoking user-level write requires open for write mode.
    check_open_for_write = function() {
      if (!private$is_open_for_write()) {
        stop(paste("Item must be open for write.", self$uri))
      }
    },

    # Per the spec, invoking user-level get-metadata requires open for
    # read mode or write mode.
    check_open_for_read_or_write = function() {
      if (!self$is_open()) {
        stop(paste("Item must be open for read or write.", self$uri))
      }
    }
  )
)
