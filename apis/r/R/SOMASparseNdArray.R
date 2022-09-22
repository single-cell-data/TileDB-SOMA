#' SOMASparseNdArray
#'
#' @description
#' `SOMASparseNdArray` is a monomorphic, sparse, N-dimensional array with offset
#' (zero-based) integer indexing on each dimension. The `SOMASparseNdArray` has
#' a user-defined schema, which includes:
#'
#' - type - a `primitive` type, expressed as an Arrow type (e.g., `int64`, `float32`, etc)
#' - shape - the shape of the array, i.e., number and length of each dimension
#'
#' All dimensions must have a positive, non-zero length.
#'
#' **Note** - on TileDB this is an sparse array with `N` uint64 dimensions of
#' domain [0, maxUint64), and a single attribute.
#'
#' ## Duplicate writes
#'
#' As duplicate index values are not allowed, index values already present in
#' the object are overwritten and new index values are added.
#'
#' @importFrom bit64 as.integer64

SOMASparseNdArray <- R6::R6Class(
  classname = "SOMASparseNdArray",
  inherit = TileDBArray,

  public = list(

    #' @description Create a SOMASparseNdArray named with the URI.
    #' @param type an [Arrow type][arrow::data-type] defining the type of each element in the array.
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
          name = paste0("__dim_", i - 1L),
          domain = bit64::as.integer64(c(0L, shape[i] - 1L)),
          tile = bit64::as.integer64(min(c(shape[i], 2048L))),
          type = "INT64"
        )
        tiledb::filter_list(tdb_dims[[i]]) <- zstd_filter_list
      }

      # create array attribute
      tdb_attr <- tiledb::tiledb_attr(
        name = "data",
        type = tiledb_type_from_arrow_type(type),
        filter_list = zstd_filter_list
      )

      # array schema
      tdb_schema <- tiledb::tiledb_array_schema(
        domain = tiledb::tiledb_domain(tdb_dims),
        attrs = tdb_attr,
        sparse = TRUE,
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

    #' @description Write matrix-like data to the array.
    #'
    #' @param data Any `matrix`-like object coercible to a
    #' [`TsparseMatrix`][`Matrix::TsparseMatrix-class`]. Character dimension
    #' names are ignored because `SOMANdArray`'s use integer indexing.
    #'
    write = function(data) {
      stopifnot(
        "'data' must be a matrix" = is_matrix(data)
      )
      # coerce to a TsparseMatrix, which uses 0-based COO indexing
      data <- as(data, Class = "TsparseMatrix")
      coo <- data.frame(
        i = bit64::as.integer64(data@i),
        j = bit64::as.integer64(data@j),
        x = data@x
      )
      colnames(coo) <- c(self$dimnames(), self$attrnames())
      private$write_coo_dataframe(coo)
    }
  ),

  private = list(

    # @description Ingest COO-formatted dataframe into the TileDB array.
    # @param x A [`data.frame`].
    write_coo_dataframe = function(data) {
      stopifnot(is.data.frame(data))
      # private$log_array_ingestion()
      on.exit(private$close())
      private$open("WRITE")
      arr <- self$object
      arr[] <- data
    }
  )
)
