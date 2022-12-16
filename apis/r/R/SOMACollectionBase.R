#' SOMA Collection Base Class
#'
#' @description Contains a key-value mapping where the keys are string names
#' and the values are any SOMA-defined foundational or composed type, including
#' `SOMACollection`, `SOMADataFrame`, `SOMADenseNdArray`, `SOMASparseNdArray`
#' or `SOMAExperiment`.
SOMACollectionBase <- R6::R6Class(
  classname = "SOMACollectionBase",
  inherit = TileDBGroup,

  public = list(

    #' @description Create a new `SOMACollection`.
    #'
    #' @param uri URI of the TileDB group
    #' @param ctx optional tiledb context
    initialize = function(uri, platform_config = NULL, ctx = NULL) {
      super$initialize(uri, platform_config, ctx)
    },

    #' @description Add a new SOMA object to the collection.
    create = function() {
      super$create()
      private$write_object_type_metadata()
    },

    #' @description Add a new SOMA object to the collection.
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

    #' @description Retrieve a SOMA object by name.
    #' @param name The name of the object to retrieve.
    #' @returns SOMA object.
    get = function(name) {
      super$get(name)
    }

  ),

  active = list(
    #' @field somas Retrieve [`SOMA`] members.
    soma_type = function(value) {
      if (!missing(value)) {
        stop("`soma_type` is a read-only field", call. = FALSE)
      }
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
      # from the object's metadata so we knowh which SOMA class to instantiate
      tiledb_constructor <- switch(type,
        ARRAY = TileDBArray$new,
        GROUP = TileDBGroup$new,
        stop(sprintf("Unknown member TileDB type: %s", type))
      )

      tiledb_object <- tiledb_constructor(uri, self$ctx, self$platform_config)
      soma_type <- tiledb_object$get_metadata(SOMA_OBJECT_TYPE_METADATA_KEY)

      soma_constructor <- switch(soma_type,
        SOMADataFrame = SOMADataFrame$new,
        SOMADenseNdArray = SOMADenseNdArray$new,
        SOMASparseNdArray = SOMASparseNdArray$new,
        SOMACollection = SOMACollection$new,
        stop(sprintf("Unknown member SOMA type: %s", soma_type))
      )
      soma_constructor(uri, self$ctx, self$platform_config)
    }
  )
)