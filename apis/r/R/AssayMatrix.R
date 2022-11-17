#' Single-cell Assay Matrix
#'
#' @description
#' Base class for 2D sparse `matrix`-like data with string dimensions. An
#' `AssayMatrix` may contain one or more "layers" (i.e., additional measurements
#' that share the same dimensions and non-empty coordinates.
#'
#' Used for the `X` field of [`SOMA`].
#' @param batch_mode logical, if `TRUE`, batch query mode is enabled, which
#' provides the ability to detect partial query results and resubmit until
#' all results are retrieved.
#' @importFrom Matrix sparseMatrix
#' @export

AssayMatrix <- R6::R6Class(
  classname = "AssayMatrix",
  inherit = AnnotationArray,

  public = list(
    #' @field verbose Print status messages
    verbose = TRUE,

    #' @description Ingest assay data from a sparse matrix
    #' @param x any `matrix`-like object coercible to a
    #' [`TsparseMatrix`][`Matrix::TsparseMatrix-class`] with string dimensions.
    #' @param index_cols Names to use for the TileDB array's dimensions that
    #' will contain the matrix row/column names.
    #' @param value_col Name to use for the TileDB array's attribute that will
    #' contain the matrix values.
    #' @param transpose If `TRUE`, the order of the new TileDB array's
    #' dimensions are reversed relative to the names defined by `index_cols`.
    #'
    #' ## Transposing Dimensions
    #'
    #' The `transpose` argument allows you to ingest a matrix into a TileDB
    #' array with the dimension order reversed relative to the original matrix
    #' shape. For example, consider matrix `mat` with dimensions `i` and `j`:
    #'
    #' ```{r}
    #' mat <- matrix(1:6, nrow = 3, dimnames = list(LETTERS[1:3], letters[1:2]))
    #' mat
    #' ```
    #'
    #' By default `AssayMatrix` will create a TileDB array will have dimensions
    #' `i`/`j` (i.e., the same order as the matrix):
    #'
    #' ```{r}
    #' amat1 <- AssayMatrix$new(uri = "mem://assaymatrix1", verbose = FALSE)
    #' amat1$from_matrix(mat, index_cols = c("i", "j"))
    #' amat1$to_dataframe()
    #' ```
    #'
    #' If you set `transpose = TRUE`, the TileDB array will have dimensions
    #' `j`/`i` instead:
    #'
    #' ```{r}
    #' amat2 <- AssayMatrix$new(uri = "mem://assaymatrix2", verbose = FALSE)
    #' amat2$from_matrix(mat, index_cols = c("i", "j"), transpose = TRUE)
    #' amat2$to_dataframe()
    #' ```

    from_matrix = function(
      x,
      index_cols,
      value_col = "value",
      transpose = FALSE
    ) {
      stopifnot(
        "Must provide 'index_cols' to name the index columns" = !missing(index_cols),
        "'value_col' must be scalar" = is_scalar_character(value_col)
      )
      private$validate_matrix(x)
      df <- matrix_to_coo(x, index_cols = index_cols, value_cols = value_col)

      # reverse index columns if transposing array dimensions
      if (transpose) index_cols <- rev(index_cols)

      self$from_dataframe(df, index_cols = index_cols)
    },

    #' @description Ingest assay data from a COO-formatted data frame
    #' @param x a [`data.frame`]
    #' @param index_cols A column index, either numeric with a column index, or
    #' character with a column name, identifying the 2 index columns. All
    #' other columns are ingested as attributes.
    from_dataframe = function(x, index_cols) {
      stopifnot(
        "Must provide 'index_cols' to identify the index columns" = !missing(index_cols),
        "'x' must be a data.frame" = is.data.frame(x),
        length(index_cols) == 2,
        all(index_cols %in% colnames(x))
      )
      if (!self$exists()) {
        private$create_empty_array(x, index_cols)
      } else {
        if (self$verbose) {
          message(
            sprintf("Updating existing %s at '%s'", self$class(), self$uri)
          )
        }
      }
      private$ingest_data(x, index_cols)
    },

    #' @description Retrieve the assay data from TileDB
    #' @param attrs Specify one or more attributes to retrieve. If `NULL`,
    #' all attributes are retrieved.
    #' @return A [`Matrix::dgTMatrix-class`].
    to_dataframe = function(attrs = NULL, batch_mode = FALSE) {
      private$read_data(
        attrs = attrs,
        batch_mode = batch_mode,
        return_as = "data.frame"
      )
    },

    #' @description Retrieve assay data from TileDB as a 2D sparse matrix.
    #' @param attr The name of the attribute layer to retrieve. If `NULL`, the
    #' first layer is returned.
    #' @param transpose If `TRUE`, the order of the matrix's dimensions are
    #' reversed relative to the TileDB array's dimension name.
    #' dimensions are reversed relative to the names defined by `index_cols`.
    #' @return A [`Matrix::dgTMatrix-class`].
    to_matrix = function(attr = NULL, batch_mode = FALSE, transpose = FALSE) {
      if (is.null(attr)) {
        attr <- self$attrnames()[1]
      }
      stopifnot(is_scalar_character(attr))

      assay_data <- private$read_data(
        attrs = attr,
        batch_mode = batch_mode,
        return_as = "data.frame"
      )

      # reverse index columns if transposing array dimensions
      if (transpose) assay_data <- assay_data[c(rev(self$dimnames()), attr)]

      assay_dims <- vapply_int(assay_data[1:2], n_unique)
      row_labels <- unique(assay_data[[1]])
      col_labels <- unique(assay_data[[2]])

      Matrix::sparseMatrix(
        i = match(assay_data[[1]], row_labels),
        j = match(assay_data[[2]], col_labels),
        x = assay_data[[3]],
        dims = assay_dims,
        dimnames = list(row_labels, col_labels),
        repr = "T"
      )
    }
  ),

  private = list(

    # @description Create an empty TileDB array with a schema optimized for 2D
    # COO-formatted data.
    create_empty_array = function(
      x,
      index_cols = c("obs_id", "var_id"),
      cell_order = "ROW_MAJOR",
      tile_order = "ROW_MAJOR",
      capacity = 1000) {

      # determine appropriate type for each attribute
      value_cols <- setdiff(colnames(x), index_cols)
      stopifnot(
        "'x' must contain >=1 non-indexing columns" = length(value_cols) >= 1
      )
      value_types <- vapply_char(x[value_cols], tiledb::r_to_tiledb_type)

      # array dimensions
      tdb_dims <- mapply(
        FUN = tiledb::tiledb_dim,
        name = index_cols,
        MoreArgs = list(
          type = "ASCII",
          domain = NULL,
          tile = NULL
        ),
        SIMPLIFY = FALSE
      )

      tiledb::filter_list(tdb_dims[[1]]) <- tiledb::tiledb_filter_list(
        tiledb::tiledb_filter("RLE")
      )

      # TODO: Make zstd compression level configurable, currently using same
      # default as core: https://github.com/TileDB-Inc/TileDB/blob/56644c1e94fcba26d07a608112fdcdf3fd120ba8/tiledb/sm/filter/compression_filter.h#L154
      tiledb::filter_list(tdb_dims[[2]]) <- tiledb::tiledb_filter_list(
        tiledb::tiledb_filter_set_option(
          object = tiledb::tiledb_filter("ZSTD"),
          option = "COMPRESSION_LEVEL",
          value = 3L
        )
      )

      # array attributes
      tdb_attr_filter <- tiledb::tiledb_filter_set_option(
        object = tiledb::tiledb_filter("ZSTD"),
        option = "COMPRESSION_LEVEL",
        value = 3L
      )

      tdb_attrs <- mapply(
        FUN = tiledb::tiledb_attr,
        name = value_cols,
        type = value_types,
        MoreArgs = list(
          filter_list = tiledb::tiledb_filter_list(tdb_attr_filter),
          ctx = self$ctx
        ),
        SIMPLIFY = FALSE
      )

      # array schema
      tdb_schema <- tiledb::tiledb_array_schema(
        domain = tiledb::tiledb_domain(tdb_dims),
        attrs = tdb_attrs,
        cell_order = cell_order,
        tile_order = tile_order,
        sparse = TRUE,
        capacity = capacity,
        offsets_filter_list = tiledb::tiledb_filter_list(c(
          tiledb::tiledb_filter("DOUBLE_DELTA"),
          tiledb::tiledb_filter("BIT_WIDTH_REDUCTION"),
          tiledb::tiledb_filter("ZSTD")
        ))
      )

      private$log_array_creation(index_cols)
      tiledb::tiledb_array_create(uri = self$uri, schema = tdb_schema)
      private$write_object_type_metadata()
    },

    # @description Ingest assay data into the TileDB array.
    # @param x A [`data.frame`] containing the assay data.
    # @param index_cols Character vector with column names to use as index
    ingest_data = function(x, index_cols) {
      stopifnot(
        "Assay data must be a data.frame" = is.data.frame(x)
      )
      private$log_array_ingestion()
      tdb_array <- tiledb::tiledb_array(self$uri, query_type = "WRITE")
      tdb_array[] <- x
    }
  )
)
