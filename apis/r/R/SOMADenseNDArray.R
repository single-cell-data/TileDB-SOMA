#' SOMA Dense Nd-Array
#'
#' @description \code{SOMADenseNDArray} is a dense, N-dimensional array of
#' a \code{primitive} type, with offset (zero-based) \code{int64} integer
#' indexing on each dimension with domain \code{[0, maxInt64)}. The
#' \code{SOMADenseNDArray} has a user-defined schema, which includes:
#' \itemize{
#'  \item \code{type}: a \code{primitive} type, expressed as an Arrow type
#'   (e.g., \code{\link[arrow]{int64}}, \code{\link[arrow]{float32}}, etc),
#'   indicating the type of data contained within the array.
#'  \item \code{shape}: the shape of the array, i.e., number and length of each
#'   dimension. This is a soft limit which can be increased using
#'   \code{$resize()} up to the \code{maxshape}.
#'  \item \code{maxshape}: the hard limit up to which \code{shape} may be
#'   increased using \code{$resize()}.
#' }
#' All dimensions must have a positive, non-zero length, and there must be 1 or
#' more dimensions.
#'
#' The default \dQuote{fill} value for \code{SOMADenseNDArray} is the zero or
#' null value of the array type (e.g.,
#' \code{\link[arrow:float32]{arrow::float32}()} defaults to 0.0).
#'
#' @param coords Optional \code{list} of integer vectors, one for each
#' dimension, with a length equal to the number of values to read. If
#' \code{NULL}, all values are read. List elements can be named when
#' specifying a subset of dimensions.
#' @template param-result-order
#' @param log_level Optional logging level with default value of
#' \dQuote{\code{warn}}.
#'
#' @export
#'
#' @inherit SOMADenseNDArrayCreate examples
#'
SOMADenseNDArray <- R6::R6Class(
  classname = "SOMADenseNDArray",
  inherit = SOMANDArrayBase,
  public = list(
    #' @description Read as an \link[arrow:Table]{Arrow table}
    #' (lifecycle: maturing).
    #'
    #' @return An \link[arrow:Table]{Arrow table}.
    #'
    read_arrow_table = function(
      coords = NULL,
      result_order = "auto",
      log_level = "auto"
    ) {
      private$.check_handle()
      stopifnot(
        "'log_level' must be character" = is.character(log_level),
        "'result_order' must be character" = is.character(result_order)
      )

      result_order <- map_query_layout(match_query_layout(result_order))

      if (is.null(coords)) {
        # These are 0-up: add 1 for R use
        ned <- self$non_empty_domain(max_only = TRUE)
        coords <- lapply(X = as.integer(ned), FUN = function(x) {
          0:x
        })
      }
      coords <- private$.convert_coords(coords)

      soma_debug(sprintf(
        "[SOMADenseNDArray$read_arrow_table] timestamp (%s)",
        self$tiledb_timestamp %||% "now"
      ))

      rl <- soma_array_read(
        soma_array = private$.handle,
        dim_points = coords,
        result_order = result_order,
        loglevel = log_level
      )

      soma_array_to_arrow_table(rl)
    },

    #' @description Read as a dense matrix (lifecycle: maturing).
    #'
    #' @return A \code{matrix}.
    #'
    read_dense_matrix = function(
      coords = NULL,
      result_order = "ROW_MAJOR",
      log_level = "warn"
    ) {
      ndim <- self$ndim()
      attrnames <- self$attrnames()

      stopifnot(
        "Array must have two dimensions" = ndim == 2,
        "Array must contain column 'soma_data'" = all.equal(
          "soma_data",
          attrnames
        )
      )

      if (is.null(coords)) {
        ned <- self$non_empty_domain(max_only = TRUE)
        # These are 0-up: add 1 for R use
        nrow <- as.numeric(ned[[1]]) + 1
        ncol <- as.numeric(ned[[2]]) + 1
      } else {
        nrow <- length(unique(as.numeric(coords[[1]])))
        ncol <- length(unique(as.numeric(coords[[2]])))
      }

      tbl <- self$read_arrow_table(
        coords = coords,
        result_order = result_order,
        log_level = log_level
      )
      return(matrix(
        as.numeric(tbl$GetColumnByName("soma_data")),
        nrow = nrow,
        ncol = ncol,
        byrow = result_order == "ROW_MAJOR"
      ))
    },

    #' @description Write matrix data to the array (lifecycle: maturing).\cr
    #' \cr
    #' \strong{Note}: The \code{$write()} method is currently limited to writing
    #' from two-dimensional matrices (lifecycle: maturing).
    #'
    # More general write methods for higher-dimensional array could be added.
    #'
    #' @param values A \code{matrix}. Character dimension names are ignored
    #' because \code{SOMANDArray}s use integer indexing.
    #' @param coords A \code{list} of integer vectors, one for each dimension,
    #' with a length equal to the number of values to write. If \code{NULL},
    #' the default, the values are taken from the row and column names
    #' of \code{values}.
    #'
    #' @return Invisibly returns \code{self}.
    #'
    write = function(values, coords = NULL) {
      private$.check_handle()
      soma_debug("[SOMADenseNDArray::write] entered")
      stopifnot("'values' must be a matrix" = is.matrix(values))

      if (is.null(coords)) {
        coords <- list(seq_len(nrow(values)), seq_len(ncol(values)))
      }

      stopifnot(
        "'coords' must be a list of integer vectors" = is.list(coords) &&
          all(vapply_lgl(coords, is.integer)),
        "length of 'coords' must match number of dimensions" = length(coords) ==
          self$ndim()
      )

      ## the 'soma_data' data type may not have been cached, and if so we need to fetch it
      if (is.null(private$.type)) {
        private$.type <- self$schema()[["soma_data"]]$type
      }

      soma_debug("[SOMADenseNDArray::write] about to call write")
      arrsch <- arrow::schema(arrow::field("soma_data", private$.type))
      tbl <- arrow::arrow_table(soma_data = values, schema = arrsch)

      soma_debug("[SOMADenseNDArray::write] array created")
      naap <- nanoarrow::nanoarrow_allocate_array()
      nasp <- nanoarrow::nanoarrow_allocate_schema()
      arrow::as_record_batch(tbl)$export_to_c(naap, nasp)
      writeArrayFromArrow(
        soma_array = private$.handle,
        naap = naap,
        nasp = nasp,
        arraytype = "SOMADenseNDArray"
      )
      soma_debug("[SOMADenseNDArray::write] written")

      return(invisible(self))
    }
  ),
  private = list(
    # @description Open the handle for the C++ interface
    .open_handle = function(open_mode, timestamp) {
      private$.handle <- open_dense_ndarray_handle(
        self$uri,
        open_mode,
        private$.context$handle,
        timestamp
      )
    },

    # Given a user-specified shape along a particular dimension, returns a named
    # list containing name, capacity, and extent elements. The shape cannot be
    # NULL for dense arrays.
    .dim_capacity_and_extent = function(name, shape = NULL, create_options) {
      out <- list(name = name, capacity = NULL, extent = NULL)
      stopifnot(
        "'shape' must be a positive scalar integer" = rlang::is_scalar_integerish(
          shape
        ) &&
          shape > 0
      )
      out$capacity <- shape
      out$extent <- min(shape, create_options$dim_tile(name))
      out
    }
  )
)
