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
      private$write_object_type_metadata()
      self
    },

    #' @description Read as an 'arrow::Table' (lifecycle: experimental)
    #' @param coords Optional `list` of integer vectors, one for each dimension, with a
    #' length equal to the number of values to read. If `NULL`, all values are
    #' read. List elements can be named when specifying a subset of dimensions.
    #' @template param-result-order
    #' @param iterated Option boolean indicated whether data is read in call (when
    #' `FALSE`, the default value) or in several iterated steps.
    #' @param log_level Optional logging level with default value of `"warn"`.
    #' @return An [`arrow::Table`].
    read_arrow_table = function(
      coords = NULL,
      result_order = "auto",
      iterated = FALSE,
      log_level = "warn"
    ) {
      uri <- self$uri

      result_order <- map_query_layout(match_query_layout(result_order))

      if (!is.null(coords)) {
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
      }

      if (isFALSE(iterated)) {
          rl <- soma_array_reader(uri = uri,
                                  dim_points = coords,        # NULL is dealt with by soma_array_reader()
                                  result_order = result_order,
                                  loglevel = log_level,       # idem
                                  config = as.character(tiledb::config(
                                      self$tiledbsoma_ctx$get_tiledb_context()
                                  )))
          private$soma_reader_transform(rl)
      } else {
          ## should we error if this isn't null?
          if (!is.null(self$soma_reader_pointer)) {
              warning("pointer not null, skipping")
          } else {
              private$soma_reader_setup()
              private$sparse_repr <- "" # no sparse matrix transformation
          }
          invisible(NULL)
      }
    },

    #' @description Read as a sparse matrix (lifecycle: experimental)
    #' @param coords Optional `list` of integer vectors, one for each dimension, with a
    #' length equal to the number of values to read. If `NULL`, all values are
    #' read. List elements can be named when specifying a subset of dimensions.
    #' @template param-result-order
    #' @param repr Optional one-character code for sparse matrix representation type
    #' @param iterated Option boolean indicated whether data is read in call (when
    #' `FALSE`, the default value) or in several iterated steps.
    #' @param log_level Optional logging level with default value of `"warn"`.
    #' @return A `matrix` object
    read_sparse_matrix = function(
      coords = NULL,
      result_order = "auto",
      repr = c("C", "T", "R"),
      iterated = FALSE,
      log_level = "warn"
    ) {
      repr <- match.arg(repr)
      dims <- self$dimensions()
      attr <- self$attributes()
      stopifnot("Array must have two dimensions" = length(dims) == 2,
                "Array must contain columns 'soma_dim_0' and 'soma_dim_1'" =
                    all.equal(c("soma_dim_0", "soma_dim_1"), names(dims)),
                "Array must contain column 'soma_data'" = all.equal("soma_data", names(attr)))

      if (isFALSE(iterated)) {
          tbl <- self$read_arrow_table(coords = coords, result_order = result_order, log_level = log_level)
          Matrix::sparseMatrix(i = 1 + as.numeric(tbl$GetColumnByName("soma_dim_0")),
                               j = 1 + as.numeric(tbl$GetColumnByName("soma_dim_1")),
                               x = as.numeric(tbl$GetColumnByName("soma_data")),
                               dims = as.integer(self$shape()), repr = repr)
      } else {
          ## should we error if this isn't null?
          if (!is.null(self$soma_reader_pointer)) {
              warning("pointer not null, skipping")
          } else {
              private$soma_reader_setup()
              private$sparse_repr <- repr
          }
          invisible(NULL)
      }
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
      stopifnot(is.data.frame(values))
      # private$log_array_ingestion()
      on.exit(self$close())
      private$open("WRITE")
      arr <- self$object
      arr[] <- values
    },

    ## refined from base class
    soma_reader_transform = function(x) {
      tbl <- as_arrow_table(x)
      if (private$sparse_repr == "") {
          tbl
      } else {
          Matrix::sparseMatrix(i = 1 + as.numeric(tbl$GetColumnByName("soma_dim_0")),
                               j = 1 + as.numeric(tbl$GetColumnByName("soma_dim_1")),
                               x = as.numeric(tbl$GetColumnByName("soma_data")),
                               dims = as.integer(self$shape()), repr = private$sparse_repr)
      }
    },

    ## internal 'repr' state variable, by default 'unset'
    sparse_repr = ""

  )
)
