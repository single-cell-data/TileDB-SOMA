#' SOMA NDArray Base Class
#'
#' @description
#' Adds NDArray-specific functionality to the [`SOMAArrayBase`] class.
#' (lifecycle: experimental)
#'
#' @export
#' @importFrom bit64 as.integer64

SOMANDArrayBase <- R6::R6Class(
  classname = "SOMANDArrayBase",
  inherit = SOMAArrayBase,

  public = list(

    #' @description Create a SOMA NDArray named with the URI. (lifecycle:
    #' experimental)
    #' @param type an [Arrow type][arrow::data-type] defining the type of each
    #' element in the array.
    #' @param shape a vector of integers defining the shape of the array.
    #' @template param-platform-config
    #' @param internal_use_only Character value to signal this is a 'permitted'
    #' call, as `create()` is considered internal and should not be called
    #' directly.
    create = function(type, shape, platform_config = NULL, internal_use_only = NULL) {
      if (is.null(internal_use_only) || internal_use_only != "allowed_use") {
        stop(paste("Use of the create() method is for internal use only. Consider using a",
                   "factory method as e.g. 'SOMASparseNDArrayCreate()'."), call. = FALSE)
      }

      tdb_schema <- private$.build_tiledb_schema(
        type = type,
        shape = shape,
        is_sparse = private$.is_sparse,
        platform_config = platform_config
      )

      # create array
      tiledb::tiledb_array_create(uri = self$uri, schema = tdb_schema)
      self$open("WRITE", internal_use_only = "allowed_use")
      private$write_object_type_metadata()
      self
    }
  ),

  private = list(
    .is_sparse = NULL,

    .build_tiledb_schema = function(
      type,
      shape,
      is_sparse,
      platform_config = NULL
    ) {

      stopifnot(
        "'type' must be a valid Arrow type" =
          is_arrow_data_type(type),
        "'shape' must be a vector of positive integers" =
          is.vector(shape) && all(shape > 0),
        "'is_sparse' must be a scalar logical" = is_scalar_logical(is_sparse)
      )

      # Parse the tiledb/create/ subkeys of the platform_config into a handy,
      # typed, queryable data structure.
      tiledb_create_options <- TileDBCreateOptions$new(platform_config)

      # Default zstd filter to use if none is specified in platform config
      default_zstd_filter <- list(
        name = "ZSTD",
        COMPRESSION_LEVEL = tiledb_create_options$dataframe_dim_zstd_level()
      )

      # create array dimensions
      tdb_dims <- vector(mode = "list", length = length(shape))
      for (i in seq_along(shape)) {
        dim_info <- private$.dim_capacity_and_extent(
          name = paste0("soma_dim_", i - 1L),
          shape = shape[i],
          create_options = tiledb_create_options
        )

        tdb_dims[[i]] <- tiledb::tiledb_dim(
          name = dim_info$name,
          domain = bit64::as.integer64(c(0L, dim_info$capacity - 1L)),
          tile = bit64::as.integer64(dim_info$extent),
          type = "INT64",
          filter_list = tiledb::tiledb_filter_list(
            filters =  tiledb_create_options$dim_filters(
              dim_name = dim_info$name,
              default = list(default_zstd_filter)
            )
          )
        )
      }

      # attribute filters
      tdb_attr_filters <- tiledb::tiledb_filter_list(
        filters = tiledb_create_options$attr_filters(
          attr_name = "soma_data",
          default = list(default_zstd_filter)
        ))

      # create array attribute
      tdb_attr <- tiledb::tiledb_attr(
        name = "soma_data",
        type = tiledb_type_from_arrow_type(type),
        filter_list = tdb_attr_filters
      )

      # array schema
      cell_tile_orders <- tiledb_create_options$cell_tile_orders()
      tiledb::tiledb_array_schema(
        domain = tiledb::tiledb_domain(tdb_dims),
        attrs = tdb_attr,
        sparse = private$.is_sparse,
        cell_order = cell_tile_orders["cell_order"],
        tile_order = cell_tile_orders["tile_order"],
        capacity = tiledb_create_options$capacity(),
        allows_dups = tiledb_create_options$allows_duplicates(),
        offsets_filter_list = tiledb::tiledb_filter_list(
          tiledb_create_options$offsets_filters()
        ),
        validity_filter_list = tiledb::tiledb_filter_list(
          tiledb_create_options$validity_filters()
        )
      )
    },

    # Given a user-specified shape along a particular dimension, returns a named
    # list containing name, capacity, and extent elements.
    .dim_capacity_and_extent = function(name, shape = NULL, create_options) {
      stop(
        ".dim_capacity_and_extent() must be implemented by child class",
        call. = FALSE
      )
    },

    #  @description Converts a list of vectors corresponding to coords to a
    #  format acceptable for sr_setup and soma_array_reader
    .convert_coords = function(coords) {

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

    # Internal marking of one or zero based matrices for iterated reads
    zero_based = NA
  )
)
