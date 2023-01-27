#' TileDB Object Base Class
#'
#' @description
#' Base class to implement shared functionality across the TileDBArray and
#' TileDBGroup classes. [lifecycle: experimental]
#' @export
TileDBObject <- R6::R6Class(
  classname = "TileDBObject",
  public = list(
    #' @field platform_config Optional platform configuration
    platform_config = NULL,
    #' @field ctx Optional TileDB context
    ctx = NULL,

    #' @description Create a new TileDB object. [lifecycle: experimental]
    #' @param uri URI for the TileDB object
    #' @param verbose Print status messages
    #' @param platform_config Optional platform configuration
    #' @param ctx optional TileDB context
    initialize = function(uri, platform_config = NULL, ctx = NULL) {
      if (missing(uri)) stop("Must specify a `uri`", call. = FALSE)
      private$tiledb_uri <- TileDBURI$new(uri)
      self$platform_config <- platform_config
      self$ctx <- ctx

      if (is.null(self$ctx)) {
        self$ctx <- tiledb::tiledb_get_context()
      }
    },

    #' @description Print the name of the R6 class. [lifecycle: experimental]
    class = function() {
      class(self)[1]
    },

    #' @description Print-friendly representation of the object. [lifecycle: experimental]
    print = function() {
      cat(glue::glue("<{self$class()}>"), sep = "\n")
      cat("  uri:", self$uri, "\n")
    },

    #' @description Check if the object exists. [lifecycle: experimental]
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
      tiledb::tiledb_object_type(self$uri, ctx = self$ctx) %in% expected_type
    }
  ),

  active = list(
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

    # @description Contains TileDBURI object [lifecycle: experimental]
    tiledb_uri = NULL

  )
)
