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
    },

    #' @param uri A URI for a SOMA object
    #' @return TRUE if the URI will use `tiledbv2` semantics.
    is_tiledbv2 = function(uri) {
      self$get_data_protocol(uri) == "tiledbv2"
    },

    #' @param uri A URI for a SOMA object
    #' @return TRUE if the URI will use `tiledbv3` semantics.
    is_tiledbv3 = function(uri) {
      self$get_data_protocol(uri) == "tiledbv3"
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

#' Set the Default Global Context
#'
#' Configure a default \code{\link{SOMAContext}} to be used by all TileDB-SOMA
#' operations when no explicit context is provided.
#'
#' This function should be called once at the beginning of your session before
#' opening any SOMA objects if you want to customize the TileDB context
#' parameters that will apply to all subsequent operations. Otherwise, a default
#' context will be created automatically with standard parameters when you first
#' open a SOMA object.
#'
#' If the global context was already set, an error will be raised unless
#' \code{replace=True}. Setting a new default context will not change the
#' context for TileDB-SOMA objects that were already created.
#'
#' @template param-config
#' @param replace Allow replacing the existing default context with new
#' configuration parameters.
#'
#' @return Invisibly, the default default context object.
#'
#' @export
#'
set_default_context <- function(config = NULL, replace = FALSE) {
  if (!replace && !is.null(.pkgenv[["somactx"]])) {
    stop(
      "A default context was already created. To replace the default context for new objects call again with `replace=True`.",
      call. = FALSE
    )
  }
  context <- SOMAContext$new(config)
  .pkgenv[["somactx"]] <- context
  invisible(context)
}

#' Get the Default SOMA Context
#' 
#' Retrieve the current default \code{\link{SOMAContext}} used by TileDB-SOMA
#' operations.
#'
#' This function returns the context that was either:
#' \itemize{
#'   \item Explicitly set via \code{\link{set_default_context}}, or
#'   \item Automatically created when a SOMA object was first created
#' }
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
    stop(
      "No default context is set. Call `set_default_context()` to initialize the context.",
      call. = FALSE
    )
  }
  return(context)
}
