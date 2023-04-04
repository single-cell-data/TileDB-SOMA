#' SOMA Array Base Class
#'
#' Adds SOMA-specific functionality to the [`TileDBArray`] class.  (lifecycle: experimental)

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

  public = list(

    #' @description Check if iterated read is complete or not. (lifecycle: experimental)
    read_complete = function() {
      if (is.null(private$soma_reader_pointer)) {
          TRUE
      } else {
          sr_complete(private$soma_reader_pointer)
      }
    },

    #' @description Read the next chunk of an iterated read. (lifecycle: experimental)
    read_next = function() {
      if (is.null(private$soma_reader_pointer)) {
          NULL
      } else {
          if (sr_complete(private$soma_reader_pointer)) {
              invisible(NULL)
          } else {
              rl <- sr_next(private$soma_reader_pointer)
              private$soma_reader_transform(rl)
          }
      }
    }

  ),

  private = list(

    # Cache object's SOMA_OBJECT_TYPE_METADATA_KEY
    soma_type_cache = NULL,

    update_soma_type_cache = function() {
      private$soma_type_cache <- self$get_metadata(SOMA_OBJECT_TYPE_METADATA_KEY)
    },

    write_object_type_metadata = function() {
      meta <- list()
      meta[[SOMA_OBJECT_TYPE_METADATA_KEY]] <- self$class()
      meta[[SOMA_ENCODING_VERSION_METADATA_KEY]] <- SOMA_ENCODING_VERSION
      self$set_metadata(meta)
    },

    # Internal 'external pointer' object used for iterated reads
    soma_reader_pointer = NULL,

    # Instantiate soma_reader_pointer with a soma_array_reader object
    soma_reader_setup = function() {
      private$soma_reader_pointer <- sr_setup(
        self$uri,
        config=as.character(tiledb::config(self$tiledbsoma_ctx$context()))
      )
    },

    ## to be refined in derived classes
    soma_reader_transform = function(x) {
      x
    }

  )
)
