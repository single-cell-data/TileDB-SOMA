#' SOMA Context
#'
#' Context map for TileDB-SOMA objects
#'
#' @export
#'
#'
SOMAContext <- R6::R6Class(
  classname = "SOMAContext",
  public = list(
    #' @template param-config
    #'
    #' @return An instantiated \code{SOMATileDBContext} object
    #'
    initialize = function(config = NULL) {
      ratio_array_data_key = "sm.mem.reader.sparse_global_order.ratio_array_data"
      if (is.null(config) || !(ratio_array_data_key %in% names(config))) {
        config[ratio_array_data_key] <- "0.3"
      }
      private$.handle <- create_soma_context(config)
    },

    #' @return A character vector with the current config options set on the context.
    #'
    get_config = function() {
      return(get_config_from_soma_context(private$.handle))
    },

    #' @param uri A URI for a SOMA object
    #'
    #' @return The data protocol to use for the URI.
    #'
    get_data_protocol = function(uri) {
      return(get_data_protocol_from_soma_context(private$.handle, uri))
    }

  ),
  active = list(
    #' @field handle External pointer to the C++ interface
    #'
    handle = function(value) {
      if (!missing(x = value)) {
        stop("Field `handle` is read-only", call. = FALSE)
      }
      return(private$.handle)
    }

  ),
  private = list(
    # @field Internal 'external pointer' for C++ interface ...
    #
    .handle = NULL
  )
)

#' Create and return a global default context object
#'
#' It is recommended to call this method once before all other TileDB-SOMA R API calls. If the global context
#' was already set, a warning will be raised. Setting a new default context will not change the context for
#' TileDB-SOMA objects that were already created.
#'
#' @template param-config
#'
#' @return The context that will be used for TileDB-SOMA API when no context is provided by the user.
#'
#' @export
#'
set_default_context <- function(config = NULL) {
  if (!is.null(.pkgenv[["somactx"]])) {
    warning("A default context was already created. Existing objects will not use the new default context.", call. = FALSE)
  }
  context <- SOMAContext$new(config)
  .pkgenv[["somactx"]] <- context
  return(context)
}

#' Returns the global default context
#'
#' An error is raised if no default context is set.
#'
#' @return The context that will be used for TileDB-SOMA API when no context is provided by the user.
#'
#' @export
#'
get_default_context <- function() {
  context <- .pkgenv[["somactx"]]
  if (is.null(context)) {
    stop("No default context is set. Call `set_default_context` to initialize the context.", call. = FALSE)
  }
  return(context)
}
