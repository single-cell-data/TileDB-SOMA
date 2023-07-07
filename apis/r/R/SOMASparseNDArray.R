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
SOMASparseNDArray <- R6::R6Class(
  classname = "SOMASparseNDArray",
  inherit = SOMANDArrayBase,

  public = list(

    #' @description Reads a user-defined slice of the \code{SOMASparseNDArray}
    #' @param coords Optional `list` of integer vectors, one for each dimension, with a
    #' length equal to the number of values to read. If `NULL`, all values are
    #' read. List elements can be named when specifying a subset of dimensions.
    #' @template param-result-order
    #' @param iterated Option boolean indicated whether data is read in call (when
    #' `FALSE`, the default value) or in several iterated steps.
    #' @param log_level Optional logging level with default value of `"warn"`.
    #' @return \link{SOMASparseNDArrayRead}
    read = function(
      coords = NULL,
      result_order = "auto",
      log_level = "auto"
    ) {
      private$check_open_for_read()
      result_order <- map_query_layout(match_query_layout(result_order))

      if (!is.null(coords)) {
        coords <- private$.convert_coords(coords)
      }

      cfg <- as.character(tiledb::config(self$tiledbsoma_ctx$context()))
      sr <- sr_setup(uri = self$uri,
                     config = cfg,
                     dim_points = coords,
                     result_order = result_order,
                     timestamp_end = private$tiledb_timestamp,
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
      private$.write_coo_dataframe(coo)
    },

    #' @description Retrieve number of non-zero elements (lifecycle: experimental)
    #' @return A scalar with the number of non-zero elements
    nnz = function() {
      nnz(self$uri, config=as.character(tiledb::config(self$tiledbsoma_ctx$context())))
    }

  ),

  private = list(
    .is_sparse = TRUE,

    # Given a user-specified shape along a particular dimension, returns a named
    # list containing name, capacity, and extent elements. If no shape is
    # provided the .Machine$integer.max - 1 is used.
    .dim_capacity_and_extent = function(name, shape = NULL, create_options) {
      out <- list(name = name, capacity = NULL, extent = NULL)

      if (is.null(shape)) {
        out$capacity <- .Machine$integer.max - 1
        out$extent <- min(out$capacity, create_options$dim_tile(name))
      } else {
        stopifnot(
          "'shape' must be a positive scalar integer" =
            rlang::is_scalar_integerish(shape) && shape > 0
        )
        out$capacity <- shape
        out$extent <- min(shape, create_options$dim_tile(name))
      }

      out
    },

    # @description Ingest COO-formatted dataframe into the TileDB array.
    # (lifecycle: experimental)
    # @param values A [`data.frame`].
    .write_coo_dataframe = function(values) {
      private$check_open_for_write()

      stopifnot(is.data.frame(values))
      # private$log_array_ingestion()
      arr <- self$object
      if (!is.null(private$tiledb_timestamp)) {
          arr@timestamp <- private$tiledb_timestamp
      }
      arr[] <- values
    },

    ## internal 'repr' state variable, by default 'unset'
    sparse_repr = "",

    # Internal marking of one or zero based matrices for iterated reads
    zero_based = NA

  )
)
