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
SOMAObject <- R6::R6Class(
  classname = "SOMAObject",
  public = list(
    #' @description Create a new SOMA object. (lifecycle: maturing)
    #'
    #' @param uri ...
    #' @param platform_config ...
    #' @param tiledbsoma_ctx ...
    #' @param tiledb_timestamp ...
    #' @param internal_use_only ...
    #' @param soma_context ...
    #'
    initialize = function(
      uri,
      platform_config = NULL,
      tiledbsoma_ctx = NULL,
      tiledb_timestamp = NULL,
      internal_use_only = NULL,
      soma_context = NULL
    ) {
      envs <- unique(vapply(
        X = unique(sys.parents()),
        FUN = function(n) environmentName(environment(sys.function(n))),
        FUN.VALUE = character(1L)
      ))
      if (!"tiledbsoma" %in% envs) {
        stop(
          paste(
            strwrap(
              paste(
                "Use of the new() method is for internal use only. Consider",
                "using a factory method such as 'SOMADataFrameOpen()'."
              ),
              width = 0.6 * getOption("width")
            ),
            collapse = '\n'
          ),
          call. = FALSE
        )
      }
      .NotYetImplemented()
    },

    #' @description Close and reopen the TileDB object in a new mode
    #'
    #' @param mode New mode to open the object in; choose from:
    #' \itemize{
    #'  \item \dQuote{\code{READ}}
    #'  \item \dQuote{\code{WRITE}}
    #' }
    #' @param tiledb_timestamp Optional Datetime (POSIXct) with TileDB timestamp
    #'
    #' @return Invisibly returns \code{self} opened in \code{mode}
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
    }

  ),
  active = list(
    #' @field platform_config Platform configuration
    #'
    platform_config = function(value) {
      if (!missing(x = value)) {
        stop("'platform_config' is a read-only field", call. = FALSE)
      }
      return(private$.tiledb_platform_config)
    },

    #' @field tiledbsoma_ctx SOMATileDBContext
    #'
    tiledbsoma_ctx = function(value) {
      if (!missing(x = value)) {
        stop("'tiledbsoma_ctx' is a read-only field", call. = FALSE)
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

    #' @field uri
    #' The URI of the TileDB object.
    uri = function(value) {
      if (!missing(value)) {
        private$.read_only_error("uri")
      }
      return(private$.uri)
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

    # @field .mode ...
    #
    .mode = character(1L),

    # @field .uri ...
    #
    .uri = character(1L),

    # @description Throw an error saying a field is read-only
    #
    .read_only_error = function(field) {
      stop("Field ", sQuote(field), " is read-only", call. = FALSE)
    }
  )
)

SOMAObjectCreate <- \() SOMAObject$new()
