#' TileDB Object Base Class
#'
#' @description
#' Base class to implement shared functionality across the TileDBArray and
#' TileDBGroup classes. (lifecycle: experimental)
#' @export
TileDBObject <- R6::R6Class(
  classname = "TileDBObject",
  public = list(
    #' @description Create a new TileDB object. (lifecycle: experimental)
    #' @param uri URI for the TileDB object
    #' @param platform_config Optional platform configuration
    #' @param tiledbsoma_ctx Optional SOMATileDBContext
    #' @param mode One of "READ" or "WRITE"
    #' @param internal_use_only Character value to signal 'permitted' call as
    #' `new()` is considered internal and should not be called directly
    initialize = function(uri, platform_config = NULL, tiledbsoma_ctx = NULL,
                          mode = "READ", internal_use_only = NULL) {
      ## calls <- vapply(
      ##   X = lapply(X = sys.calls(), FUN = as.character),
      ##   FUN = '[[',
      ##   FUN.VALUE = character(length = 1L),
      ##   1L
      ## )
      ## if ('TileDBObject$new' %in% calls) {
      ##   .NotYetImplemented()
      ## }
      if (is.null(internal_use_only) || internal_use_only != "allowed_use") {
        stop(paste("Use of the new() method is discouraged, consider using a",
                   "factory method as e.g. 'SOMADataFrameOpen()'."), call. = FALSE)
      }
      if (missing(uri)) stop("Must specify a `uri`", call. = FALSE)
      private$tiledb_uri <- TileDBURI$new(uri)
      # Set platform config
      platform_config <- platform_config %||% PlatformConfig$new()
      if (!inherits(platform_config, 'PlatformConfig')) {
        stop("'platform_config' must be a PlatformConfig object", call. = FALSE)
      }
      private$tiledb_platform_config <- platform_config
      # Set context
      tiledbsoma_ctx <- tiledbsoma_ctx %||% SOMATileDBContext$new()
      if (!inherits(x = tiledbsoma_ctx, what = 'SOMATileDBContext')) {
        stop("'tiledbsoma_ctx' must be a SOMATileDBContext object", call. = FALSE)
      }
      private$.tiledbsoma_ctx <- tiledbsoma_ctx

      private$mode <- match.arg(mode, c("READ", "WRITE"))
      spdl::debug("[TileDBObject] initialize {} with '{}' in {}", self$class(), self$uri, private$mode)
    },

    #' @description Print the name of the R6 class.
    class = function() {
      class(self)[1]
    },

    #' @description Print-friendly representation of the object.
    print = function() {
      cat(glue::glue("<{self$class()}>"), sep = "\n")
      cat("  uri:", self$uri, "\n")
    },

    #' @description Check if the object exists. (lifecycle: experimental)
    #' @return `TRUE`` if the object exists, `FALSE` otherwise.
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
      tiledb::tiledb_object_type(self$uri, ctx = self$tiledbsoma_ctx$context()) %in% expected_type
    },
    #' @param param Parameter name from \code{self$platform_config} to fetch
    #' @return SOMATileDBContext
    get_tiledb_config = function(param = NULL) {
      if (!is.null(x = param)) {
        cfg <- suppressWarnings(expr = self$platform_config$get(
          platform = 'tiledb',
          param = param,
          default = NULL
        ))
        if (inherits(x = cfg, what = 'MappingBase')) {
          private$.tiledbsoma_ctx$update(cfg)
        }
      }
      return(self$tiledbsoma_ctx)
    }
  ),

  active = list(
    #' @field platform_config Platform configuration
    platform_config = function(value) {
      if (!missing(x = value)) {
        stop("'platform_config' is a read-only field", call. = FALSE)
      }
      return(private$tiledb_platform_config)
    },
    #' @field tiledbsoma_ctx SOMATileDBContext
    tiledbsoma_ctx = function(value) {
      if (!missing(x = value)) {
        stop("'tiledbsoma_ctx' is a read-only field", call. = FALSE)
      }
      return(private$.tiledbsoma_ctx)
    },
    #' @field uri
    #' The URI of the TileDB object.
    uri = function(value) {
      if (missing(value)) return(private$tiledb_uri$uri)
      stop(sprintf("'%s' is a read-only field.", "uri"), call. = FALSE)
    },

    #' @field object Access the underlying TileB object directly (either a
    #' [`tiledb::tiledb_array`] or [`tiledb::tiledb_group`]).
    object = function(value) {
      if (!missing(value)) {
        stop(sprintf("'%s' is a read-only field.", "object"), call. = FALSE)
      }
      # If the array was created after the object was instantiated, we need to
      # initialize private$tiledb_object
      if (is.null(private$tiledb_object)) {
        private$initialize_object()
      }
      private$tiledb_object
    }
  ),

  private = list(

    # Internal pointer to the TileDB object
    tiledb_object = NULL,

    # @description Contains TileDBURI object
    tiledb_uri = NULL,

    # Internal platform config
    tiledb_platform_config = NULL,

    # Internal context
    .tiledbsoma_ctx = NULL,

    # Internal mode: one of READ or WRITE
    mode = NULL,

    .read_only_error = function(field) {
      stop("Field ", sQuote(field), " is read-only", call. = FALSE)
    }

  )
)
