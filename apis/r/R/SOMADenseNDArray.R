#' SOMADenseNDArray
#'
#' @description
#' `SOMADenseNDArray` is a dense, N-dimensional array of `primitive` type, with
#' offset (zero-based) `int64` integer indexing on each dimension with domain
#' `[0, maxInt64)`. The `SOMADenseNDArray` has a user-defined schema, which
#' includes:
#'
#' - **type**: a `primitive` type, expressed as an Arrow type (e.g., `int64`,
#'   `float32`, etc), indicating the type of data contained within the array
#' - **shape**: the shape of the array, i.e., number and length of each
#'   dimension
#'
#' All dimensions must have a positive, non-zero length, and there must be 1 or
#' more dimensions.
#'
#' The default "fill" value for `SOMADenseNDArray` is the zero or null value of
#' the array type (e.g., Arrow.float32 defaults to 0.0).
#'
#' The `write` method is currently limited to writing from 2-d matrices.
#' (lifecycle: experimental)
#' @export

SOMADenseNDArray <- R6::R6Class(
  classname = "SOMADenseNDArray",
  inherit = SOMAArrayBase,

  public = list(

    #' @description Create a SOMADenseNDArray named with the URI. (lifecycle: experimental)
    #' @param type an [Arrow type][arrow::data-type] defining the type of each
    #' element in the array.
    #' @param shape a vector of integers defining the shape of the array.
    #' @template param-platform-config
    #' @param internal_use_only Character value to signal this is a 'permitted' call,
    #' as `create()` is considered internal and should not be called directly.
    create = function(type, shape, platform_config=NULL, internal_use_only = NULL) {
      if (is.null(internal_use_only) || internal_use_only != "allowed_use") {
        stop(paste("Use of the create() method is for internal use only. Consider using a",
                   "factory method as e.g. 'SOMADenseNDArrayCreate()'."), call. = FALSE)
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
      # use tiledb default names like `__dim_0`
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
            list(
              list(name="ZSTD", COMPRESSION_LEVEL=tiledb_create_options$dataframe_dim_zstd_level())
            )
          )
        )
      }

      # create array attribute
      tdb_attr <- tiledb::tiledb_attr(
        name = "soma_data",
        type = tiledb_type_from_arrow_type(type),
        filter_list = tiledb::tiledb_filter_list(tiledb_create_options$attr_filters("soma_data"))
      )

      # array schema
      cell_tile_orders <- tiledb_create_options$cell_tile_orders()
      tdb_schema <- tiledb::tiledb_array_schema(
        domain = tiledb::tiledb_domain(tdb_dims),
        attrs = tdb_attr,
        sparse = FALSE,
        cell_order = cell_tile_orders["cell_order"],
        tile_order = cell_tile_orders["tile_order"],
        capacity=tiledb_create_options$capacity(),
        allows_dups=tiledb_create_options$allows_duplicates(),
        offsets_filter_list = tiledb::tiledb_filter_list(
          tiledb_create_options$offsets_filters()
        ),
        validity_filter_list = tiledb::tiledb_filter_list(
          tiledb_create_options$validity_filters()
        )
      )

      # create array
      tiledb::tiledb_array_create(uri = self$uri, schema = tdb_schema)
      self$open("WRITE", internal_use_only = "allowed_use")
      private$write_object_type_metadata()
      self
    },

    #' @description Read as an 'arrow::Table' (lifecycle: experimental)
    #' @param coords Optional `list` of integer vectors, one for each dimension, with a
    #' length equal to the number of values to read. If `NULL`, all values are
    #' read. List elements can be named when specifying a subset of dimensions.
    #' @template param-result-order
    #' @param log_level Optional logging level with default value of `"warn"`.
    #' @return An [`arrow::Table`].
    read_arrow_table = function(
      coords = NULL,
      result_order = "auto",
      log_level = "warn"
    ) {
      private$check_open_for_read()

      uri <- self$uri

      result_order <- map_query_layout(match_query_layout(result_order))

      if (!is.null(coords)) {
        coords <- private$convert_coords(coords)
      }

      cfg <- as.character(tiledb::config(self$tiledbsoma_ctx$context()))
      rl <- soma_array_reader(uri = uri,
                              dim_points = coords,        # NULL dealt with by soma_array_reader()
                              result_order = result_order,
                              loglevel = log_level,       # idem
                              config = cfg)

      soma_array_to_arrow_table(rl)
    },

    #' @description Read as a dense matrix (lifecycle: experimental)
    #' @param coords Optional `list` of integer vectors, one for each dimension, with a
    #' length equal to the number of values to read. If `NULL`, all values are
    #' read. List elements can be named when specifying a subset of dimensions.
    #' @template param-result-order
    #' @param log_level Optional logging level with default value of `"warn"`.
    #' @return A `matrix` object
    read_dense_matrix = function(
      coords = NULL,
      result_order = "ROW_MAJOR",
      log_level = "warn"
    ) {
      private$check_open_for_read()

      dims <- self$dimensions()
      attr <- self$attributes()
      stopifnot("Array must have two dimensions" = length(dims) == 2,
                "Array must contain columns 'soma_dim_0' and 'soma_dim_1'" =
                    all.equal(c("soma_dim_0", "soma_dim_1"), names(dims)),
                "Array must contain column 'soma_data'" = all.equal("soma_data", names(attr)))

      tbl <- self$read_arrow_table(coords = coords, result_order = result_order, log_level = log_level)
      m <- matrix(as.numeric(tbl$GetColumnByName("soma_data")),
                  nrow = length(unique(as.numeric(tbl$GetColumnByName("soma_dim_0")))),
                  ncol = length(unique(as.numeric(tbl$GetColumnByName("soma_dim_1")))),
                  byrow = result_order == "ROW_MAJOR")

    },

    #' @description Write matrix data to the array. (lifecycle: experimental)
    #'
    #' More general write methods for higher-dimensional array could be added.
    #'
    #' @param values A `matrix`. Character dimension names are ignored because
    #' `SOMANDArray`'s use integer indexing.
    #' @param coords A `list` of integer vectors, one for each dimension, with a
    #' length equal to the number of values to write. If `NULL`, the default,
    #' the values are taken from the row and column names of `values`.
    write = function(values, coords = NULL) {
      private$check_open_for_write()

      stopifnot(
        "'values' must be a matrix" = is.matrix(values)
      )

      if (is.null(coords)) {
        coords <- list(seq_len(nrow(values)), seq_len(ncol(values)))
      }

      stopifnot(
        "'coords' must be a list of integer vectors" =
          is.list(coords) && all(vapply_lgl(coords, is.integer)),
        "length of 'coords' must match number of dimensions" =
          length(coords) == length(self$dimensions())
      )

      arr <- self$object
      tiledb::query_layout(arr) <- "COL_MAJOR"
      arr[] <- values
    }
  ),

  private = list(

    #  @description Converts a list of vectors corresponding to coords to a
    #  format acceptable for sr_setup and soma_array_reader
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
