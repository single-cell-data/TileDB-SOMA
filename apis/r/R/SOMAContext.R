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
