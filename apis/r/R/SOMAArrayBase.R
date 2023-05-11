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
    
    #' @description Converts a list of vectors corresponding to coords to a 
    #' format acceptable for sr_setup and soma_array_reader
    convert_coords = function(coords) {
      
      ## ensure coords is a named list, use to select dim points
      stopifnot("'coords' must be a list" = is.list(coords),
                "'coords' must be a list of vectors or integer64" =
                    all(vapply_lgl(coords, is_vector_or_int64)),
                "'coords' if unnamed must have length of dim names, else if named names must match dim names" =
                    (is.null(names(coords)) && length(coords) == length(self$dimnames())) ||
                    (!is.null(names(coords)) && all(names(coords) %in% self$dimnames()))
                )

      ## if unnamed (and test for length has passed in previous statement) set names
      if (is.null(names(coords))) names(coords) <- self$dimnames()

      ## convert integer to integer64 to match dimension type
      coords <- lapply(coords, function(x) if (inherits(x, "integer")) bit64::as.integer64(x) else x)
      
      coords
      
    }

  )
)
