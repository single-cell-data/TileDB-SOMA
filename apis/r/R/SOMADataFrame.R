#' SOMADataFrame
#'
#' @description
#' `SOMADataFrame` is a multi-column table containing a "pseudo-column" called
#' `soma_rowid`, of type `uint64` and domain `[0, #rows)`. `soma_rowid` is the
#' row offset (row id), and is a contiguous integer number beginning with zero.
#'
#' `SOMADataFrame` must contain a column called `soma_joinid`, of type `uint64`.
#' The `soma_joinid` column contains a unique value for each row in the
#' `SOMADataFrame`, and intended to act as a joint key for other objects, such
#' as [`SOMASparseNdArray`].
#' @importFrom stats setNames
#' @export

SOMADataFrame <- R6::R6Class(
  classname = "SOMADataFrame",
  inherit = TileDBArray,

  public = list(

    #' @description Create
    #' @param schema an [`arrow::schema`].
    create = function(schema) {
      stopifnot(
        "'schema' must be a valid Arrow schema" =
          inherits(schema, "Schema"),
        "'soma_rowid' is a reserved column name" =
          !"soma_rowid" %in% schema$names
      )

      # soma_rowid array dimension
      tdb_dim <- tiledb::tiledb_dim(
        name = "soma_rowid",
        # TODO: tiledbsoma-py uses the full uint64 range but we can't do that
        # here (see https://issues.apache.org/jira/browse/ARROW-9083)
        domain = bit64::as.integer64(c(0, 2^32 - 1)),
        tile = bit64::as.integer64(2048),
        type = "UINT64",
        filter_list = tiledb::tiledb_filter_list(c(
          tiledb_zstd_filter()
        ))
      )

      # create array attributes
      tdb_attrs <- stats::setNames(
        object = vector(mode = "list", length = schema$num_fields),
        nm = schema$names
      )
      for (field_name in schema$names) {
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
        domain = tiledb::tiledb_domain(tdb_dim),
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
    #' @param values An [`arrow::RecordBatch`] containing all columns, including
    #' the `soma_rowid` index column. The schema for `values` must match the
    #' schema for the `SOMADataFrame`.
    #'
    write = function(values) {
      stopifnot(
        "'values' must be an Arrow RecordBatch" =
          is_arrow_record_batch(values),
        "'values' must contain a 'soma_rowid' column name" =
          "soma_rowid" %in% values$names()
      )

      df <- as.data.frame(values)[c(self$dimnames(), self$attrnames())]
      # coerce from integer to integer64 in order to match the dimension type
      df$soma_rowid <- bit64::as.integer64(df$soma_rowid)

      on.exit(private$close())
      private$open("WRITE")
      arr <- self$object
      arr[] <- df
    }
  )
)
