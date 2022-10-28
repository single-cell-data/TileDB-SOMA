#' SOMADenseNdArray
#'
#' @description
#' `SOMADenseNdArray` is a dense, N-dimensional array of `primitive` type, with
#' offset (zero-based) `int64` integer indexing on each dimension with domain
#' `[0, maxUint64)`. The `SOMADenseNdArray` has a user-defined schema, which
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
#' The default "fill" value for `SOMADenseNdArray` is the zero or null value of
#' the array type (e.g., Arrow.float32 defaults to 0.0).
SOMADenseNdArray <- R6::R6Class(
  classname = "SOMADenseNdArray",
  inherit = TileDBArray,

  public = list(

    #' @description Create a SOMADenseNdArray named with the URI.
    #' @param type an [Arrow type][arrow::data-type] defining the type of each
    #' element in the array.
    #' @param shape a vector of integers defining the shape of the array.
    create = function(type, shape) {
      stopifnot(
        "'type' must be a valid Arrow type" =
          is_arrow_data_type(type),
        "'shape' must be a vector of positive integers" =
          is.vector(shape) && all(shape > 0)
      )

      zstd_filter_list <- tiledb::tiledb_filter_list(c(
          tiledb_zstd_filter(level = 3L)
      ))

      # create array dimensions
      # use tiledb default names like `__dim_0`
      tdb_dims <- vector(mode = "list", length = length(shape))
      for (i in seq_along(shape)) {
        tdb_dims[[i]] <- tiledb::tiledb_dim(
          name = paste0("soma_dim_", i - 1L),
          domain = bit64::as.integer64(c(0L, shape[i] - 1L)),
          tile = bit64::as.integer64(min(c(shape[i], 2048L))),
          type = "INT64"
        )
        tiledb::filter_list(tdb_dims[[i]]) <- zstd_filter_list
      }

      # create array attribute
      tdb_attr <- tiledb::tiledb_attr(
        name = "soma_data",
        type = tiledb_type_from_arrow_type(type),
        filter_list = zstd_filter_list
      )

      # array schema
      tdb_schema <- tiledb::tiledb_array_schema(
        domain = tiledb::tiledb_domain(tdb_dims),
        attrs = tdb_attr,
        sparse = FALSE,
        cell_order = "ROW_MAJOR",
        tile_order = "ROW_MAJOR",
        capacity=100000,
        offsets_filter_list = tiledb::tiledb_filter_list(c(
          tiledb::tiledb_filter("DOUBLE_DELTA"),
          tiledb::tiledb_filter("BIT_WIDTH_REDUCTION"),
          tiledb::tiledb_filter("ZSTD")
        ))
      )

      # create array
      tiledb::tiledb_array_create(uri = self$uri, schema = tdb_schema)
    },

    #' @description Read as an ['arrow::Table']
    #' @param coords A `list` of integer vectors, one for each dimension, with a
    #' length equal to the number of values to read. If `NULL`, all values are
    #' read.
    #' @param result_order Order of read results. This can be one of either
    #' `"ROW_MAJOR, `"COL_MAJOR"`, `"GLOBAL_ORDER"`, or `"UNORDERED"`.
    #' @return An [`arrow::Table`].
    read_arrow_table = function(
      coords = NULL,
      result_order = "ROW_MAJOR"
    ) {
      on.exit(private$close())
      private$open("READ")

      arr <- self$object

      # select ranges
      if (!is.null(coords)) {
        tiledb::selected_ranges(arr) <- coords
      }

      # result order
      tiledb::query_layout(arr) <- match_query_layout(result_order)
      tiledb::return_as(arr) <- "asis"

      do.call(arrow::arrow_table, arr[])
    },

    #' @description Write matrix data to the array.
    #'
    #' @param values A `matrix`. Character dimension names are ignored because
    #' `SOMANdArray`'s use integer indexing.
    #' @param coords A `list` of integer vectors, one for each dimension, with a
    #' length equal to the number of values to write. If `NULL`, the default,
    #' the values are taken from the row and column names of `values`.
    write = function(values, coords = NULL) {
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

      on.exit(private$close())
      private$open("WRITE")
      arr <- self$object
      tiledb::query_layout(arr) <- "COL_MAJOR"
      arr[] <- values
    }
  )
)
