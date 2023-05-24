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
      private$check_open_for_read()

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
          cfg <- as.character(tiledb::config(self$tiledbsoma_ctx$context()))
          rl <- soma_array_reader(uri = uri,
                                  dim_points = coords,        # NULL dealt with by soma_array_reader()
                                  result_order = result_order,
                                  loglevel = log_level,       # idem
                                  config = cfg)
          private$soma_reader_transform(rl)
      } else {
          ## should we error if this isn't null?
          if (!is.null(self$soma_reader_pointer)) {
              warning("Reader pointer not null, skipping")
              rl <- NULL
          } else {
              private$soma_reader_setup()
              private$sparse_repr <- "" # no sparse matrix transformation
              rl <- list()
              while (!self$read_complete()) {
                  ## soma_reader_transform() applied inside read_next()
                  rl <- c(rl, self$read_next())
              }
          }
          invisible(rl)
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
    #' @return A \link{matrixZeroBasedView}
    read_sparse_matrix_zero_based = function(
      coords = NULL,
      result_order = "auto",
      repr = c("C", "T", "R"),
      iterated = FALSE,
      log_level = "warn"
    ) {
      private$check_open_for_read()

      repr <- match.arg(repr)
      dims <- self$dimensions()
      attr <- self$attributes()
      stopifnot("Array must have two dimensions" = length(dims) == 2,
                "Array must contain columns 'soma_dim_0' and 'soma_dim_1'" =
                    all.equal(c("soma_dim_0", "soma_dim_1"), names(dims)),
                "Array must contain column 'soma_data'" = all.equal("soma_data", names(attr)))

      if (isFALSE(iterated)) {
          tbl <- self$read_arrow_table(coords = coords, result_order = result_order, log_level = log_level)
          # To instantiate the one-based Matrix::sparseMatrix, we need to add 1 to the
          # zero-based soma_dim_0 and soma_dim_1. But, because these dimensions are
          # usually populated with soma_joinid, users will need to access the matrix
          # using the original, possibly-zero IDs. Therefore, we'll wrap the one-based
          # sparseMatrix with a shim providing basic access with zero-based indexes.
          # If needed, user can then explicitly ask the shim for the underlying
          # sparseMatrix using `as.one.based()`.
          mat <- Matrix::sparseMatrix(i = 1 + as.numeric(tbl$GetColumnByName("soma_dim_0")),
                                      j = 1 + as.numeric(tbl$GetColumnByName("soma_dim_1")),
                                      x = as.numeric(tbl$GetColumnByName("soma_data")),
                                      dims = as.integer(self$shape()), repr = repr)
          matrixZeroBasedView$new(mat)
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
      private$check_open_for_write()

      stopifnot(is.data.frame(values))
      # private$log_array_ingestion()
      arr <- self$object
      arr[] <- values
    },

    ## refined from base class
    soma_reader_transform = function(x) {
      tbl <- as_arrow_table(x)
      if (private$sparse_repr == "") {
          tbl
      } else {
          mat <- Matrix::sparseMatrix(i = 1 + as.numeric(tbl$GetColumnByName("soma_dim_0")),
                                      j = 1 + as.numeric(tbl$GetColumnByName("soma_dim_1")),
                                      x = as.numeric(tbl$GetColumnByName("soma_data")),
                                      dims = as.integer(self$shape()), repr = private$sparse_repr)
          # see read_sparse_matrix_zero_based() abave
          matrixZeroBasedView$new(mat)
      }
    },

    ## internal 'repr' state variable, by default 'unset'
    sparse_repr = ""

  )
)
