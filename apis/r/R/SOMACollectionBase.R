#' SOMA Collection Base Class
#'
#' @description Base class for objects containing persistent collection of SOMA
#' objects, mapping string keys to any SOMA object.  (lifecycle: experimental)
SOMACollectionBase <- R6::R6Class(
  classname = "SOMACollectionBase",
  inherit = TileDBGroup,

  public = list(

    #' @description Create a new `SOMACollection`. (lifecycle: experimental)
    #'
    #' @param uri URI of the TileDB group
    #' @param platform_config Optional storage-engine specific configuration
    #' @param ctx optional tiledb context
    initialize = function(uri, platform_config = NULL, ctx = NULL) {
      super$initialize(uri, platform_config, ctx)
    },

    #' @description Add a new SOMA object to the collection. (lifecycle: experimental)
    create = function() {
      super$create()
      private$write_object_type_metadata()
      self
    },

    #' @description Add a new SOMA object to the collection. (lifecycle: experimental)
    #' @param object SOMA object.
    #' @param name The name to use for the object. Defaults to the object URI's
    #' base name.
    #' @param relative An optional logical value indicating whether the new
    #' object's URI is relative to the collection's URI. If `NULL` (the
    #' default), the object's URI is assumed to be relative unless it is a
    #' `tiledb://` URI.
    set = function(object, name = NULL, relative = NULL) {
      # TODO: Check that object is a SOMA object
      super$set(object, name, relative)
    },

    #' @description Retrieve a SOMA object by name. (lifecycle: experimental)
    #' @param name The name of the object to retrieve.
    #' @returns SOMA object.
    get = function(name) {
      super$get(name)
    }

  ),

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

    # Instantiate a soma member object.
    # Responsible for calling the appropriate R6 class constructor.
    construct_member = function(uri, type) {
      stopifnot(
        is_scalar_character(uri),
        is_scalar_character(type)
      )

      # We have to use the appropriate TileDB base class to read the soma_type
      # from the object's metadata so we know which SOMA class to instantiate
      tiledb_constructor <- switch(type,
        ARRAY = TileDBArray$new,
        GROUP = TileDBGroup$new,
        stop(sprintf("Unknown member TileDB type: %s", type), call. = FALSE)
      )

      tiledb_object <- tiledb_constructor(uri, self$ctx, self$platform_config)
      soma_type <- tiledb_object$get_metadata(SOMA_OBJECT_TYPE_METADATA_KEY)

      soma_constructor <- switch(soma_type,
        SOMADataFrame = SOMADataFrame$new,
        SOMADenseNDArray = SOMADenseNDArray$new,
        SOMASparseNDArray = SOMASparseNDArray$new,
        SOMACollection = SOMACollection$new,
        SOMAMeasurement = SOMAMeasurement$new,
        SOMAExperiment = SOMAExperiment$new,
        stop(sprintf("Unknown member SOMA type: %s", soma_type), call. = FALSE)
      )
      soma_constructor(uri, self$ctx, self$platform_config)
    },

    # Internal method called by SOMA Measurement/Experiment's active bindings
    # to retrieve or set one of the pre-defined SOMA fields (e.g., obs, X, etc). (lifecycle: experimental)
    # @param value the optional argument passed to the active binding.
    # @param name the name of the field to retrieve or set.
    # @param expected_class the expected class of the value to set.
    get_or_set_soma_field = function(value, name, expected_class) {
      if (missing(value)) return(self$get(name))

      stopifnot(
        "Must define 'name' of the field to set" = !missing(name),
        "Must define the field's 'expected_class'" = !missing(expected_class)
      )

      if (!inherits(value, expected_class)) {
        stop(
          sprintf("%s must be a '%s'", name, expected_class),
          call. = FALSE
        )
      }
      self$set(value, name = name)
    }
  )
)
