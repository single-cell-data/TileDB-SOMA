#' SOMASparseNDArray
#'
#' @description
#' `SOMASparseNDArray` is a sparse, N-dimensional array with offset
#' (zero-based) integer indexing on each dimension. The `SOMASparseNDArray` has
#' a user-defined schema, which includes:
#'
#' - type - a `primitive` type, expressed as an Arrow type (e.g., `int64`, `float32`, etc)
#' - shape - the shape of the array, i.e., number and length of each dimension
#'
#' All dimensions must have a positive, non-zero length.
#'
#' **Note** - on TileDB this is an sparse array with `N` int64 dimensions of
#' domain [0, maxInt64), and a single attribute.
#'
#' ## Duplicate writes
#'
#' As duplicate index values are not allowed, index values already present in
#' the object are overwritten and new index values are added. (lifecycle: experimental)
#'
#' @export
#' @importFrom bit64 as.integer64

SOMASparseNDArray <- R6::R6Class(
  classname = "SOMASparseNDArray",
  inherit = SOMAArrayBase,

  public = list(

    #' @description Create a SOMASparseNDArray named with the URI. (lifecycle: experimental)
    #' @param type an [Arrow type][arrow::data-type] defining the type of each element in the array.
    #' @param shape a vector of integers defining the shape of the array.
    #' @template param-platform-config
    #' @param internal_use_only Character value to signal this is a 'permitted' call,
    #' as `create()` is considered internal and should not be called directly.
    create = function(type, shape, platform_config=NULL, internal_use_only = NULL) {
      if (is.null(internal_use_only) || internal_use_only != "allowed_use") {
        stop(paste("Use of the create() method is for internal use only. Consider using a",
                   "factory method as e.g. 'SOMASparseNDArrayCreate()'."), call. = FALSE)
      }

      stopifnot(
        "'type' must be a valid Arrow type" =
          is_arrow_data_type(type),
        "'shape' must be a vector of positive integers" =
          is.vector(shape) && all(shape > 0)
      )

      # Parse the tiledb/create/ subkeys of the platform_config into a handy,
      # typed, queryable data structure.
      tiledb_create_options <- TileDBCreateOptions$new(platform_config)

      # create array dimensions
      tdb_dims <- vector(mode = "list", length = length(shape))
      for (i in seq_along(shape)) {
        dim_name <- paste0("soma_dim_", i - 1L)
        tile_extent <- tiledb_create_options$dim_tile(dim_name)
        tile_extent <- bit64::as.integer64(min(c(shape[i], tile_extent)))

        tdb_dims[[i]] <- tiledb::tiledb_dim(
          name = dim_name,
          domain = bit64::as.integer64(c(0L, shape[i] - 1L)),
          tile = tile_extent,
          type = "INT64"
        )
        tiledb::filter_list(tdb_dims[[i]]) <- tiledb::tiledb_filter_list(
          tiledb_create_options$dim_filters(
            dim_name,
            # Default to use if there is nothing specified in tiledb-create options
            # in the platform config:
            list(list(name="ZSTD", COMPRESSION_LEVEL=tiledb_create_options$dataframe_dim_zstd_level()))
          )
        )
      }

      # create array attribute
      tdb_attr <- tiledb::tiledb_attr(
        name = "soma_data",
        type = tiledb_type_from_arrow_type(type),
        filter_list = tiledb::tiledb_filter_list(tiledb_create_options$attr_filters(
          "soma_data",
          # Default to use if there is nothing specified in tiledb-create options
          # in the platform config:
          list(list(name="ZSTD", COMPRESSION_LEVEL=tiledb_create_options$dataframe_dim_zstd_level()))
        ))
      )

      # array schema
      cell_tile_orders <- tiledb_create_options$cell_tile_orders()
      tdb_schema <- tiledb::tiledb_array_schema(
        domain = tiledb::tiledb_domain(tdb_dims),
        attrs = tdb_attr,
        sparse = TRUE,
        cell_order = cell_tile_orders["cell_order"],
        tile_order = cell_tile_orders["tile_order"],
        capacity=tiledb_create_options$capacity(),
        allows_dups=tiledb_create_options$allows_duplicates(),
        offsets_filter_list = tiledb::tiledb_filter_list(tiledb_create_options$offsets_filters()),
        validity_filter_list = tiledb::tiledb_filter_list(tiledb_create_options$validity_filters())
      )

      # create array
      tiledb::tiledb_array_create(uri = self$uri, schema = tdb_schema)
      self$open("WRITE", internal_use_only = "allowed_use")
      private$write_object_type_metadata()
      self
    },

    #' @description Reads a user-defined slice of the \code{SOMASparseNDArray}
    #' @param coords Optional `list` of integer vectors, one for each dimension, with a
    #' length equal to the number of values to read. If `NULL`, all values are
    #' read. List elements can be named when specifying a subset of dimensions.
    #' @template param-result-order
    #' @param iterated Option boolean indicated whether data is read in call (when
    #' `FALSE`, the default value) or in several iterated steps.
    #' @param log_level Optional logging level with default value of `"warn"`.
    #' @return arrow::\link[arrow]{Table} or \link{TableReadIter}
    read = function(
      coords = NULL,
      result_order = "auto",
      log_level = "warn"
    ) {
      private$check_open_for_read()

      uri <- self$uri
      
      if (self$nnz() > .Machine$integer.max) {
          warning("Iteration results cannot be concatenated on its entirerity beceause ",
                  "array has non-zero elements greater than '.Machine$integer.max'.")
      }

      result_order <- map_query_layout(match_query_layout(result_order))

      if (!is.null(coords)) {
        coords <- private$convert_coords(coords)
      }

      cfg <- as.character(tiledb::config(self$tiledbsoma_ctx$context()))
      sr <- sr_setup(uri = uri, 
                     config = cfg, 
                     dim_points = coords, 
                     #result_order = result_order,
                     loglevel = log_level)
      
      SOMASparseNDArrayRead$new(sr, shape = self$shape())
    },

    #' @description Write matrix-like data to the array. (lifecycle: experimental)
    #'
    #' @param values Any `matrix`-like object coercible to a
    #' [`TsparseMatrix`][`Matrix::TsparseMatrix-class`]. Character dimension
    #' names are ignored because `SOMANDArray`'s use integer indexing.
    #'
    write = function(values) {
      stopifnot(
        "'values' must be a matrix" = is_matrix(values)
      )
      # coerce to a TsparseMatrix, which uses 0-based COO indexing
      values <- as(values, Class = "TsparseMatrix")
      coo <- data.frame(
        i = bit64::as.integer64(values@i),
        j = bit64::as.integer64(values@j),
        x = values@x
      )
      colnames(coo) <- c(self$dimnames(), self$attrnames())
      private$write_coo_dataframe(coo)
    },

    #' @description Retrieve number of non-zero elements (lifecycle: experimental)
    #' @return A scalar with the number of non-zero elements
    nnz = function() {
      nnz(self$uri, config=as.character(tiledb::config(self$tiledbsoma_ctx$context())))
    }

  ),

  private = list(

    # @description Ingest COO-formatted dataframe into the TileDB array. (lifecycle: experimental)
    # @param x A [`data.frame`].
    write_coo_dataframe = function(values) {
      private$check_open_for_write()

      stopifnot(is.data.frame(values))
      # private$log_array_ingestion()
      arr <- self$object
      arr[] <- values
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
    },

    ## internal 'repr' state variable, by default 'unset'
    sparse_repr = "",

    # Internal marking of one or zero based matrices for iterated reads
    zero_based = NA

  )
)
