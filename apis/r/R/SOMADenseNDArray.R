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
#' (lifecycle: maturing)
#' @export
SOMADenseNDArray <- R6::R6Class(
  classname = "SOMADenseNDArray",
  inherit = SOMANDArrayBase,

  public = list(

    #' @description Read as an 'arrow::Table' (lifecycle: maturing)
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

      if (is.null(coords)) {
        # These are 0-up: add 1 for R use
        ned <- self$non_empty_domain()
        coords <- lapply(X=as.integer(ned), FUN=function(x){0:x})
      }
      coords <- private$.convert_coords(coords)

      soma_context <- soma_context()  # package-level cached instance
      spdl::debug(
        "[SOMADenseNDArray$read_arrow_table] timestamp ({})",
        self$tiledb_timestamp %||% "now"
      )

      rl <- soma_array_reader(uri = uri,
                              dim_points = coords,
                              result_order = result_order,
                              timestamprange = self$.tiledb_timestamp_range,
                              soma_context = soma_context,
                              loglevel = log_level)

      soma_array_to_arrow_table(rl)
    },

    #' @description Read as a dense matrix (lifecycle: maturing)
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
                "Array must contain column 'soma_data'" = all.equal("soma_data", names(attr)))

      if (is.null(coords)) {
        ned <- self$non_empty_domain()
        # These are 0-up: add 1 for R use
        nrow <- as.numeric(ned[[1]]) + 1
        ncol <- as.numeric(ned[[2]]) + 1
      } else {
        nrow <- length(unique(as.numeric(coords[[1]])))
        ncol <- length(unique(as.numeric(coords[[2]])))
      }

      tbl <- self$read_arrow_table(coords = coords, result_order = result_order, log_level = log_level)
      m <- matrix(as.numeric(tbl$GetColumnByName("soma_data")),
                  nrow = nrow, ncol = ncol,
                  byrow = result_order == "ROW_MAJOR")

    },

    #' @description Write matrix data to the array. (lifecycle: maturing)
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

      spdl::debug("[SOMADenseNDArray::write] entered")
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

      ## the 'soma_data' data type may not have been cached, and if so we need to fetch it
      if (is.null(private$.type)) {
          ## TODO: replace with a libtiledbsoma accessor as discussed
          tpstr <- tiledb::datatype(tiledb::attrs(tiledb::schema(self$uri))[["soma_data"]])
          arstr <- arrow_type_from_tiledb_type(tpstr)
          private$.type <- arstr
      }

      arr <- self$object
      tiledb::query_layout(arr) <- "COL_MAJOR"
      spdl::debug("[SOMADenseNDArray::write] about to call write")
      arrsch <- arrow::schema(arrow::field("soma_data", private$.type))
      tbl <- arrow::arrow_table(soma_data = values, schema = arrsch)

      spdl::debug("[SOMADenseNDArray::write] array created")
      naap <- nanoarrow::nanoarrow_allocate_array()
      nasp <- nanoarrow::nanoarrow_allocate_schema()
      arrow::as_record_batch(tbl)$export_to_c(naap, nasp)
      #arr[] <- values
      writeArrayFromArrow(
        uri = self$uri,
        naap = naap,
        nasp = nasp,
        ctxxp = private$.soma_context,
        arraytype = "SOMADenseNDArray",
        config = NULL,
        tsvec = self$.tiledb_timestamp_range
      )
      spdl::debug("[SOMADenseNDArray::write] written")

      # tiledb-r always closes the array after a write operation so we need to
      # manually reopen it until close-on-write is optional
      #self$open("WRITE", internal_use_only = "allowed_use")
      invisible(self)
    }
  ),

  private = list(
    .is_sparse = FALSE,

    # Given a user-specified shape along a particular dimension, returns a named
    # list containing name, capacity, and extent elements. The shape cannot be
    # NULL for dense arrays.
    .dim_capacity_and_extent = function(name, shape = NULL, create_options) {
      out <- list(name = name, capacity = NULL, extent = NULL)
      stopifnot(
        "'shape' must be a positive scalar integer" =
          rlang::is_scalar_integerish(shape) && shape > 0
      )
      out$capacity <- shape
      out$extent <- min(shape, create_options$dim_tile(name))
      out
    }
  )
)
