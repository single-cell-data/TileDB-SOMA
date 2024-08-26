#' SOMADataFrame
#'
#' @description
#' `SOMADataFrame` is a multi-column table that must contain a column
#' called `soma_joinid` of type `int64`, which contains a unique value for each
#' row and is intended to act as a join key for other objects, such as
#' [`SOMASparseNDArray`].  (lifecycle: maturing)

#' @importFrom stats setNames
#' @export

SOMADataFrame <- R6::R6Class(
  classname = "SOMADataFrame",
  inherit = SOMAArrayBase,

  public = list(

    #' @description Create (lifecycle: maturing)
    #' @param schema an [`arrow::schema`].
    #' @param index_column_names A vector of column names to use as user-defined
    #' index columns.  All named columns must exist in the schema, and at least
    #' one index column name is required.
    #' @template param-platform-config
    #' @param timestamps Optional timestamp start and end range
    #' @param internal_use_only Character value to signal this is a 'permitted' call,
    #' as `create()` is considered internal and should not be called directly.
    create = function(schema, index_column_names = c("soma_joinid"),
                      platform_config = NULL, timestamps = NULL, internal_use_only = NULL) {
      if (is.null(internal_use_only) || internal_use_only != "allowed_use") {
        stop(paste("Use of the create() method is for internal use only. Consider using a",
                   "factory method as e.g. 'SOMADataFrameCreate()'."), call. = FALSE)
      }
      schema <- private$validate_schema(schema, index_column_names)

      attr_column_names <- setdiff(schema$names, index_column_names)
      stopifnot("At least one non-index column must be defined in the schema" =
                length(attr_column_names) > 0)

      # Parse the tiledb/create/ subkeys of the platform_config into a handy,
      # typed, queryable data structure.
      tiledb_create_options <- TileDBCreateOptions$new(platform_config)

      ## we (currently pass domain and extent values in an arrow table (i.e. data.frame alike)
      ## where each dimension is one column (of the same type as in the schema) followed by three
      ## values for the domain pair and the extent
      dom_ext_tbl <- get_domain_and_extent_dataframe(schema, index_column_names, tiledb_create_options)

      ## we transfer to the arrow table via a pair of array and schema pointers
      dnaap <- nanoarrow::nanoarrow_allocate_array()
      dnasp <- nanoarrow::nanoarrow_allocate_schema()
      arrow::as_record_batch(dom_ext_tbl)$export_to_c(dnaap, dnasp)

      ## we need a schema pointer to transfer the schema information
      nasp <- nanoarrow::nanoarrow_allocate_schema()
      schema$export_to_c(nasp)

      spdl::debug("[SOMADataFrame$create] about to create schema from arrow")
      ctxptr <- super$tiledbsoma_ctx$context()
      createSchemaFromArrow(uri = self$uri, nasp, dnaap, dnasp, TRUE, "SOMADataFrame",
                            tiledb_create_options$to_list(FALSE), soma_context(), timestamps)

      spdl::debug("[SOMADataFrame$create] about to call write_object_type_metadata")
      private$write_object_type_metadata(timestamps)

      self$open("WRITE", internal_use_only = "allowed_use")
      self
    },

    #' @description Write (lifecycle: maturing)
    #'
    #' @param values An [`arrow::Table`] or [`arrow::RecordBatch`]
    #' containing all columns, including any index columns. The
    #' schema for `values` must match the schema for the `SOMADataFrame`.
    #'
    write = function(values) {
      private$check_open_for_write()

      # Prevent downcasting of int64 to int32 when materializing a column
      op <- options(arrow.int64_downcast = FALSE)
      on.exit(options(op), add = TRUE, after = FALSE)

      schema_names <- c(self$dimnames(), self$attrnames())
      col_names <- if (is_arrow_record_batch(values)) {
                       arrow::as_arrow_table(values)$ColumnNames()
                   } else {
                       values$ColumnNames()
                   }
      stopifnot(
        "'values' must be an Arrow Table or RecordBatch" =
          (is_arrow_table(values) || is_arrow_record_batch(values)),
        "All columns in 'values' must be defined in the schema" =
          all(col_names %in% schema_names),
        "All schema fields must be present in 'values'" =
          all(schema_names %in% col_names)
      )

      ## we transfer to the arrow table via a pair of array and schema pointers
      naap <- nanoarrow::nanoarrow_allocate_array()
      nasp <- nanoarrow::nanoarrow_allocate_schema()
      arrow::as_record_batch(values)$export_to_c(naap, nasp)

      df <- as.data.frame(values)[schema_names]
      arr <- self$object
      writeArrayFromArrow(
        uri = self$uri,
        naap = naap,
        nasp = nasp,
        arraytype = "SOMADataFrame",
        config = NULL,
        tsvec = self$.tiledb_timestamp_range
      )

      invisible(self)
    },

    #' @description Read (lifecycle: maturing)
    #' Read a user-defined subset of data, addressed by the dataframe indexing
    #' column, and optionally filtered.
    #' @param coords Optional named list of indices specifying the rows to read; each (named)
    #' list element corresponds to a dimension of the same name.
    #' @param column_names Optional character vector of column names to return.
    #' @param value_filter Optional string containing a logical expression that is used
    #' to filter the returned values. See [`tiledb::parse_query_condition`] for
    #' more information.
    #' @template param-result-order
    #' @param iterated Option boolean indicated whether data is read in call (when
    #' `FALSE`, the default value) or in several iterated steps.
    #' @param log_level Optional logging level with default value of `"warn"`.
    #' @return arrow::\link[arrow]{Table} or \link{TableReadIter}
    read = function(coords = NULL,
                    column_names = NULL,
                    value_filter = NULL,
                    result_order = "auto",
                    iterated = FALSE,
                    log_level = "auto") {

      private$check_open_for_read()

      result_order <- match_query_layout(result_order)
      uri <- self$uri
      arr <- self$object                 # need array (schema) to properly parse query condition

      ## if unnamed set names
      if (!is.null(coords)) {
          if (!is.list(coords))
              coords <- list(coords)
          if (is.null(names(coords)))
              names(coords) <- self$dimnames()
      }

      stopifnot(
          ## check columns
          "'column_names' must only contain valid dimension or attribute columns" =
              is.null(column_names) || all(column_names %in% c(self$dimnames(), self$attrnames()))
      )

      coords <- validate_read_coords(coords, dimnames = self$dimnames(), schema = self$schema())

      if (!is.null(value_filter)) {
          value_filter <- validate_read_value_filter(value_filter)
          parsed <- do.call(what = tiledb::parse_query_condition,
                            args = list(expr = str2lang(value_filter), ta = arr))
          value_filter <- parsed@ptr
      }
      spdl::debug("[SOMADataFrame$read] calling sr_setup for {} at ({},{})", self$uri,
                  private$tiledb_timestamp[1], private$tiledb_timestamp[2])
      cfg <- as.character(tiledb::config(self$tiledbsoma_ctx$context()))
      rl <- sr_setup(uri = self$uri,
                     config = cfg,
                     colnames = column_names,
                     qc = value_filter,
                     dim_points = coords,
                     timestamprange = self$.tiledb_timestamp_range,  # NULL or two-elem vector
                     loglevel = log_level)
      private$ctx_ptr <- rl$ctx
      TableReadIter$new(rl$sr)
    },

    #' @description Update (lifecycle: maturing)
    #' @details
    #' Update the existing `SOMADataFrame` to add or remove columns based on the
    #' input:
    #' - columns present in the current the `SOMADataFrame` but absent from the
    #'   new `values` will be dropped
    #' - columns absent in current `SOMADataFrame` but present in the new
    #'   `values` will be added
    #' - any columns present in both will be left alone, with the exception that
    #'   if `values` has a different type for the column, the entire update
    #'   will fail because attribute types cannot be changed.
    #'
    #' Furthermore, `values` must contain the same number of rows as the current
    #' `SOMADataFrame`.
    #'
    #' @param values A `data.frame`, [`arrow::Table`], or
    #' [`arrow::RecordBatch`].
    #' @param row_index_name An optional scalar character. If provided, and if
    #' the `values` argument is a `data.frame` with row names, then the row
    #' names will be extracted and added as a new column to the `data.frame`
    #' prior to performing the update. The name of this new column will be set
    #' to the value specified by `row_index_name`.
    update = function(values, row_index_name = NULL) {
      private$check_open_for_write()
      stopifnot(
        "'values' must be a data.frame, Arrow Table or RecordBatch" =
          is.data.frame(values) || is_arrow_table(values) || is_arrow_record_batch(values)
      )

      if (is.data.frame(values)) {
        if (!is.null(row_index_name)) {
          stopifnot(
            "'row_index_name' must be NULL or a scalar character vector" =
              is_scalar_character(row_index_name),
            "'row_index_name' conflicts with an existing column name" =
              !row_index_name %in% colnames(values)
          )
          values[[row_index_name]] <- rownames(values)
        }
        values <- arrow::as_arrow_table(values)
      }

      # Retrieve existing soma_joinids from array to:
      # - validate number of rows in values matches number of rows in array
      # - add original soma_joinids to values if not present
      spdl::debug("[SOMADataFrame update]: Retrieving existing soma_joinids")
      self$reopen(mode = "READ")
      joinids <- self$read(column_names = "soma_joinid")$concat()$soma_joinid
      if (length(joinids) != nrow(values)) {
        stop(
          "Number of rows in 'values' must match number of rows in array:\n",
          "  - Number of rows in array: ", length(joinids), "\n",
          "  - Number of rows in 'values': ", nrow(values),
          call. = FALSE
        )
      }

      # Add soma_joinid column if not present
      if (!"soma_joinid" %in% colnames(values)) {
        values$soma_joinid <- joinids
      }
      private$validate_schema(
        schema = values$schema,
        index_column_names = self$dimnames()
      )

      old_schema <- self$schema()
      new_schema <- values$schema

      old_cols <- old_schema$names
      new_cols <- new_schema$names

      drop_cols <- setdiff(old_cols, new_cols)
      add_cols <- setdiff(new_cols, old_cols)
      common_cols <- intersect(old_cols, new_cols)

      tiledb_create_options <- TileDBCreateOptions$new(self$platform_config)

      # Check compatibility of new/old data types in common columns
      check_arrow_schema_data_types(
        old_schema[common_cols],
        new_schema[common_cols]
      )

      # Drop columns
      se <- tiledb::tiledb_array_schema_evolution()
      for (drop_col in drop_cols) {
        spdl::info("[SOMADataFrame update]: dropping column '{}'", drop_col)
        se <- tiledb::tiledb_array_schema_evolution_drop_attribute(
          object = se,
          attrname = drop_col
        )
      }

      # Add columns
      for (add_col in add_cols) {
        spdl::info("[SOMADataFrame update]: adding column '{}'", add_col)

        col_type <- new_schema$GetFieldByName(add_col)$type
        attr <- tiledb_attr_from_arrow_field(
          field = new_schema$GetFieldByName(add_col),
          tiledb_create_options = tiledb_create_options
        )

        if (inherits(col_type, "DictionaryType")) {
          spdl::debug(
            "[SOMADataFrame update]: adding column '{}' as an enumerated type",
            add_col
          )
          se <- tiledb::tiledb_array_schema_evolution_add_enumeration(
            object = se,
            name = add_col,
            enums = levels(values$GetColumnByName(add_col)$as_vector()),
            ordered = col_type$ordered
          )
          attr <- tiledb::tiledb_attribute_set_enumeration_name(attr, add_col)
        }

        se <- tiledb::tiledb_array_schema_evolution_add_attribute(se, attr)
      }

      se <- tiledb::tiledb_array_schema_evolution_array_evolve(se, self$uri)

      # Reopen array for writing with new schema
      self$reopen(mode = "WRITE")
      spdl::info("[SOMADataFrame update]: Writing new data")
      self$write(values)
    },

    #' @description Retrieve the shape; as \code{SOMADataFrames} are shapeless,
    #' simply raises an error
    #'
    #' @return None, instead a \code{\link{.NotYetImplemented}()} error is raised
    #'
    shape = function() stop(errorCondition(
      "'SOMADataFrame$shape()' is not implemented yet",
      class = 'notYetImplementedError'
    ))

  ),

  private = list(

    # @description Validate schema (lifecycle: maturing)
    # Handle default column additions (eg, soma_joinid) and error checking on
    # required columns
    # @return An [`arrow::Schema`], which may be modified by the addition of
    # required columns.
    validate_schema = function(schema, index_column_names) {
      stopifnot(
        "'schema' must be a valid Arrow schema" =
          is_arrow_schema(schema),
        "'index_column_names' must be a non-empty character vector" =
            is.character(index_column_names) && length(index_column_names) > 0,
        "All 'index_column_names' must be defined in the 'schema'" =
          assert_subset(index_column_names, schema$names, "indexed field"),
        "Column names must not start with reserved prefix 'soma_'" =
          all(!startsWith(setdiff(schema$names, "soma_joinid"), "soma_"))
      )

      # Add soma_joinid column if not present
      if ("soma_joinid" %in% schema$names) {
        stopifnot(
          "soma_joinid field must be of type Arrow int64" =
            schema$GetFieldByName("soma_joinid")$type == arrow::int64()
        )
      } else {
        schema <- schema$AddField(
          i = 0,
          field = arrow::field("soma_joinid", arrow::int64())
        )
      }

      schema
    },

    # Internal variable to hold onto context returned by sr_setup
    ctx_ptr = NULL
  )
)
