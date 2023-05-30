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
    #' @param tiledbsoma_ctx optional SOMATileDBContext
    #' @param tiledb_timestamp Optional Datetime (POSIXct) for TileDB timestamp
    #' @param internal_use_only Character value to signal this is a 'permitted' call,
    #' as `new()` is considered internal and should not be called directly.
    initialize = function(uri, platform_config = NULL, tiledbsoma_ctx = NULL, tiledb_timestamp = NULL,
                          internal_use_only = NULL) {
      super$initialize(uri=uri, platform_config=platform_config,
                       tiledbsoma_ctx=tiledbsoma_ctx, tiledb_timestamp = tiledb_timestamp,
                       internal_use_only=internal_use_only)
    },

    #' @description Add a new SOMA object to the collection. (lifecycle: experimental)
    #' @param internal_use_only Character value to signal this is a 'permitted' call,
    #' as `create()` is considered internal and should not be called directly.
    create = function(internal_use_only = NULL) {
      if (is.null(internal_use_only) || internal_use_only != "allowed_use") {
        stop(paste("Use of the create() method is for internal use only. Consider using a",
                   "factory method as e.g. 'SOMACollectionCreate()'."), call. = FALSE)
      }
      super$create(internal_use_only=internal_use_only)

      # Root SOMA objects include a `dataset_type` entry to allow the
      # TileDB Cloud UI to detect that they are SOMA datasets.
      if (self$class() == "SOMAExperiment") {
        metadata <- list(dataset_type = "soma")
      } else {
        metadata <- list()
      }

      # Key constraint on the TileDB Group API:
      # * tiledb_group_create is synchronous: new items appear on disk once that
      #   function has returned.
      # * When we have an open handle for write and we do tiledb_group_put_metadata,
      #   metadata items do _not_ appear on disk until tiledb_group_close has been
      #   called.
      # To avoid confusion, we pay the price of a close-then-reopen here.
      #
      # TODO: this feels like it could be resolved through improved caching.
      self$open("WRITE", internal_use_only = "allowed_use")
      private$write_object_type_metadata(metadata)
      self$close()
      self$open(mode = "WRITE", internal_use_only = "allowed_use")
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
    #' @param mode Mode to open in
    #' @returns SOMA object.
    get = function(name) {
      super$get(name)
    },

    #' @description Add a new SOMA collection to this collection. (lifecycle: experimental)
    #' @param object SOMA collection object.
    #' @param key The key to be added.
    add_new_collection = function(object, key) {
      # TODO: Check that object is a collection
      super$set(object, key)
      object
    },

    #' @description Add a new SOMA dataframe to this collection. (lifecycle: experimental)
    #' @param key The key to be added.
    #' @param schema Arrow schema argument passed on to DataFrame$create()
    #' @param index_column_names Index column names passed on to DataFrame$create()
    add_new_dataframe = function(key, schema, index_column_names) {
      ## TODO: Check argument validity
      ## TODO: platform_config ?
      ndf <- SOMADataFrame$new(file.path(self$uri, key), internal_use_only = "allowed_use")
      ndf$create(schema, index_column_names, internal_use_only = "allowed_use")
      super$set(ndf, key)
      ndf
    },

    #' @description Add a new SOMA DenseNdArray to this collection. (lifecycle: experimental)
    #' @param key The key to be added.
    #' @param type an [Arrow type][arrow::data-type] defining the type of each
    #' element in the array.
    #' @param shape a vector of integers defining the shape of the array.
    add_new_dense_ndarray = function(key, type, shape) {
      ndarr <- SOMADenseNDArray$new(file.path(self$uri, key), internal_use_only = "allowed_use")
      ndarr$create(type, shape, internal_use_only = "allowed_use")
      super$set(ndarr, key)
      ndarr
    },

    #' @description Add a new SOMA SparseNdArray to this collection. (lifecycle: experimental)
    #' @param key The key to be added.
    #' @param type an [Arrow type][arrow::data-type] defining the type of each
    #' element in the array.
    #' @param shape a vector of integers defining the shape of the array.
    add_new_sparse_ndarray = function(key, type, shape) {
      ndarr <- SOMASparseNDArray$new(file.path(self$uri, key), internal_use_only = "allowed_use")
      ndarr$create(type, shape, internal_use_only = "allowed_use")
      super$set(ndarr, key)
      ndarr
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

    # Add standard SOMA metadata values to the object with the option to include
    # additional metadata.
    write_object_type_metadata = function(metadata = list()) {
      stopifnot(is.list(metadata))
      metadata[[SOMA_OBJECT_TYPE_METADATA_KEY]] <- self$class()
      metadata[[SOMA_ENCODING_VERSION_METADATA_KEY]] <- SOMA_ENCODING_VERSION
      self$set_metadata(metadata)
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
      tiledbsoma_constructor <- switch(type,
        ARRAY = TileDBArray$new,
        GROUP = TileDBGroup$new,
        stop(sprintf("Unknown member TileDB type: %s", type), call. = FALSE)
      )

      tiledb_object <- tiledbsoma_constructor(
        uri,
        tiledbsoma_ctx = self$tiledbsoma_ctx,
        platform_config = self$platform_config,
        tiledb_timestamp = self$.group_open_timestamp,
        internal_use_only = "allowed_use"
      )

      tiledb_object$open(mode = "READ", internal_use_only = "allowed_use")
      soma_type <- tiledb_object$get_metadata(SOMA_OBJECT_TYPE_METADATA_KEY)
      tiledb_object$close()

      spdl::debug(
        "[SOMACollectionBase] Instantiating {} object at: '{}'",
        soma_type %||% "Unknown",
        uri
      )

      stopifnot("Discovered metadata object type is missing; cannot construct" = !is.null(soma_type))
      soma_constructor <- switch(soma_type,
        SOMADataFrame = SOMADataFrame$new,
        SOMADenseNDArray = SOMADenseNDArray$new,
        SOMASparseNDArray = SOMASparseNDArray$new,
        SOMACollection = SOMACollection$new,
        SOMAMeasurement = SOMAMeasurement$new,
        SOMAExperiment = SOMAExperiment$new,
        stop(sprintf("Unknown member SOMA type: %s", soma_type), call. = FALSE)
      )
      obj <- soma_constructor(
        uri,
        tiledbsoma_ctx = self$tiledbsoma_ctx,
        platform_config = self$platform_config,
        tiledb_timestamp = self$.group_open_timestamp,
        internal_use_only = "allowed_use"
      )
      obj
    },

    # Internal method called by SOMA Measurement/Experiment's active bindings
    # to retrieve or set one of the pre-defined SOMA fields (e.g., obs, X, etc). (lifecycle: experimental)
    # @param value the optional argument passed to the active binding.
    # @param name the name of the field to retrieve or set.
    # @param expected_class the expected class of the value to set.
    get_or_set_soma_field = function(value, name, expected_class) {
      private$check_open_for_read_or_write()

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
