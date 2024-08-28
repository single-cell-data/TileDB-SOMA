#' SOMA Array Base Class
#'
#' Adds SOMA-specific functionality to the [`TileDBArray`] class.  (lifecycle:
#' experimental)
#' @keywords internal

SOMAArrayBase <- R6::R6Class(
  classname = "SOMAArrayBase",
  inherit = TileDBArray,

  active = list(
    #' @field soma_type Retrieve the SOMA object type.
    soma_type = function(value) {
      stopifnot("'soma_type' is a read-only field" = missing(value))
      if (is.null(private$soma_type_cache)) {
        private$update_soma_type_cache()
      }
      private$soma_type_cache
    }
  ),


  private = list(

    # Cache object's SOMA_OBJECT_TYPE_METADATA_KEY
    soma_type_cache = NULL,

    update_soma_type_cache = function() {
      private$soma_type_cache <- self$get_metadata(SOMA_OBJECT_TYPE_METADATA_KEY)
    },

    write_object_type_metadata = function() {
      #private$check_open_for_write()

      meta <- list()
      meta[[SOMA_OBJECT_TYPE_METADATA_KEY]] <- self$class()
      meta[[SOMA_ENCODING_VERSION_METADATA_KEY]] <- SOMA_ENCODING_VERSION
      spdl::debug("[SOMAArrayBase::write_object_metadata] calling set metadata")
      self$set_metadata(meta)
    }
  )
)
