#' SOMADataFrame
#'
#' @description
#' `SOMAIndexedDataFrame` is a multi-column table that must contain a column
#' called `soma_joinid` of type `uint64`, which contains a unique value for each
#' row and is intended to act as a join key for other objects, such as
#' [`SOMASparseNdArray`].

#' @importFrom stats setNames
#' @export

SOMAIndexedDataFrame <- R6::R6Class(
  classname = "SOMAIndexedDataFrame",
  inherit = TileDBArray,

  public = list(

    #' @description Create
    #' @param schema an [`arrow::schema`].
    create = function(schema, index_column_names) {
      stopifnot(
        "'schema' must be a valid Arrow schema" =
          inherits(schema, "Schema"),
        "'soma_rowid' is a reserved column name" =
          !"soma_rowid" %in% schema$names,
        is.character(index_column_names) && length(index_column_names) > 0,
        "All 'index_column_names' must be defined in the 'schema'" =
          (index_column_names %in% schema$names)
      )

      attr_column_names <- setdiff(schema$names, index_column_names)
      stopifnot(
        "At least one non-index column must be defined in the schema" =
          length(attr_column_names) > 0
      )

      # array dimensions
      tdb_dims <- stats::setNames(
        object = vector(mode = "list", length = length(index_column_names)),
        nm = index_column_names
      )

      for (field_name in index_column_names) {
        field <- schema$GetFieldByName(field_name)

        tdb_dims[[field_name]] <- tiledb::tiledb_dim(
          name = field_name,
          domain = arrow_type_range(field$type),
          # TODO: Parameterize
          tile = 2048L,
          type = tiledb_type_from_arrow_type(field$type),
          filter_list = tiledb::tiledb_filter_list(c(
            tiledb_zstd_filter()
          ))
        )
      }

      # array attributes
      tdb_attrs <- stats::setNames(
        object = vector(mode = "list", length = length(attr_column_names)),
        nm = attr_column_names
      )

      for (field_name in attr_column_names) {
        field <- schema$GetFieldByName(field_name)
        field_type <- tiledb_type_from_arrow_type(field$type)

        tdb_attrs[[field_name]] <- tiledb::tiledb_attr(
          name = field_name,
          type = field_type,
          nullable = field$nullable,
          ncells = if (field_type == "ASCII") NA_integer_ else 1L,
          filter_list = tiledb::tiledb_filter_list(c(
            tiledb_zstd_filter()
          ))
        )
      }

      # array schema
      tdb_schema <- tiledb::tiledb_array_schema(
        domain = tiledb::tiledb_domain(tdb_dims),
        attrs = tdb_attrs,
        sparse = TRUE,
        cell_order = "ROW_MAJOR",
        tile_order = "ROW_MAJOR",
        capacity = 100000,
        # TODO: should be configurable via a global option
        allows_dups = FALSE,
        offsets_filter_list = tiledb::tiledb_filter_list(c(
          tiledb::tiledb_filter("DOUBLE_DELTA"),
          tiledb::tiledb_filter("BIT_WIDTH_REDUCTION"),
          tiledb_zstd_filter()
        ))
      )

      # create array
      tiledb::tiledb_array_create(uri = self$uri, schema = tdb_schema)
    },

    #' @description Write
    #'
    #' @param values An [`arrow::Table`] containing all columns, including
    #' the `soma_rowid` index column. The schema for `values` must match the
    #' schema for the `SOMAIndexedDataFrame`.
    #'
    write = function(values) {
      on.exit(private$close())
      schema_names <- c(self$dimnames(), self$attrnames())

      stopifnot(
        "'values' must be an Arrow Table" =
          is_arrow_table(values),
        "All columns in 'values' must be defined in the schema" =
          all(values$ColumnNames() %in% schema_names),
        "All schema fields must be present in 'values'" =
          all(schema_names %in% values$ColumnNames())
      )

      df <- as.data.frame(values)[schema_names]
      private$open("WRITE")
      arr <- self$object
      arr[] <- df
    },

    #' @description Read
    #' Read a user-defined subset of data, addressed by the dataframe indexing
    #' column, and optionally filtered.
    #' @param ids Indices specifying the rows to read.
    #' @param column_names Character vector of column names to return.
    #' @param value_filter A string containing a logical expression that is used
    #' to filter the returned values. See [`tiledb::parse_query_condition`] for
    #' more information.
    #' @param result_order Order of read results. This can be one of either
    #' `"ROW_MAJOR, `"COL_MAJOR"`, `"GLOBAL_ORDER"`, or `"UNORDERED"`.
    #' @return An [`arrow::Table`].
    read = function(
      ids = NULL,
      column_names = NULL,
      value_filter = NULL,
      result_order = "UNORDERED"
    ) {
      on.exit(private$close())
      private$open("READ")

      arr <- self$object

      # select columns
      if (!is.null(column_names)) {
        stopifnot(
          "'column_names' must only contain non-index columns" =
            all(!column_names %in% self$dimnames())
        )
        tiledb::attrs(arr) <- column_names
      }

      # select ranges
      if (!is.null(ids)) {
        tiledb::selected_ranges(arr) <- list(cbind(ids, ids))
      }

      # filter
      if (!is.null(value_filter)) {
        stopifnot(is_scalar_character(value_filter))
        tiledb::query_condition(arr) <- do.call(
          what = tiledb::parse_query_condition,
          args = list(expr = str2lang(value_filter), ta = self$object)
        )
      }

      # result order
      tiledb::query_layout(arr) <- match_query_layout(result_order)

      arrow::arrow_table(as.data.frame(arr[]))
    }
  )
)
