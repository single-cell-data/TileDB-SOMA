#' SOMADataFrame
#'
#' @description
#' `SOMADataFrame` is a multi-column table that must contain a column
#' called `soma_joinid` of type `int64`, which contains a unique value for each
#' row and is intended to act as a join key for other objects, such as
#' [`SOMASparseNDArray`].  (lifecycle: experimental)

#' @importFrom stats setNames
#' @export

SOMADataFrame <- R6::R6Class(
  classname = "SOMADataFrame",
  inherit = SOMAArrayBase,

  public = list(

    #' @description Create (lifecycle: experimental)
    #' @param schema an [`arrow::schema`].
    #' @param index_column_names A vector of column names to use as user-defined
    #' index columns.  All named columns must exist in the schema, and at least
    #' one index column name is required.
    #' @template param-platform-config
    #' @param internal_use_only Character value to signal this is a 'permitted' call,
    #' as `create()` is considered internal and should not be called directly.
    create = function(schema, index_column_names = c("soma_joinid"),
                        platform_config = NULL, internal_use_only = NULL) {
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

      # array dimensions
      tdb_dims <- stats::setNames(
        object = vector(mode = "list", length = length(index_column_names)),
        nm = index_column_names
      )

      for (field_name in index_column_names) {
        field <- schema$GetFieldByName(field_name)

        tile_extent <- tiledb_create_options$dim_tile(field_name)

        tile_extent <- switch(field$type$ToString(),
          "int8" = as.integer(tile_extent),
          "int16" = as.integer(tile_extent),
          "int32" = as.integer(tile_extent),
          "int64" = bit64::as.integer64(tile_extent),
          "double" = as.double(tile_extent),
          "string" = NULL,
          tile_extent
        )

        # Default 2048 mods to 0 for 8-bit types and 0 is an invalid extent
        if (field$type$bit_width %||% 0L == 8L) {
          tile_extent <- 64L
        }

        tdb_opts <- tiledb_create_options$dim_filters(field_name,
          ## Default if there is nothing specified in tiledb-create options in the platform config:
          list(list(name="ZSTD", COMPRESSION_LEVEL=tiledb_create_options$dataframe_dim_zstd_level()))
        )
        tdb_dims[[field_name]] <- tiledb::tiledb_dim(
          name = field_name,
          # Numeric index types must be positive values for indexing
          domain = arrow_type_unsigned_range(field$type), tile = tile_extent,
          type = tiledb_type_from_arrow_type(field$type, is_dim=TRUE),
          filter_list = tiledb::tiledb_filter_list(tdb_opts)
        )
      }

      # array attributes
      tdb_attrs <- stats::setNames(
        object = vector(mode = "list", length = length(attr_column_names)),
        nm = attr_column_names
      )

      for (field_name in attr_column_names) {
        field <- schema$GetFieldByName(field_name)
        field_type <- tiledb_type_from_arrow_type(field$type, is_dim=FALSE)

        ## # Check if the field is ordered and mark it as such
        ## if (!is.null(x = levels[[field_name]]) && isTRUE(field$type$ordered)) {
        ##   attr(levels[[field_name]], 'ordered') <- attr(levels[[field_name]], 'ordered', exact = TRUE) %||% TRUE
        ## }

        tdb_attrs[[field_name]] <- tiledb::tiledb_attr(
          name = field_name,
          type = field_type,
          nullable = field$nullable,
          ncells = if (field_type == "ASCII" || field_type == "UTF8" ) NA_integer_ else 1L,
          filter_list = tiledb::tiledb_filter_list(tiledb_create_options$attr_filters(field_name))
        )
      }

      # array schema
      cell_tile_orders <- tiledb_create_options$cell_tile_orders()
      tdb_schema <- tiledb::tiledb_array_schema(
        domain = tiledb::tiledb_domain(tdb_dims),
        attrs = tdb_attrs,
        sparse = TRUE,
        cell_order = cell_tile_orders["cell_order"],
        tile_order = cell_tile_orders["tile_order"],
        capacity = tiledb_create_options$capacity(),
        allows_dups = tiledb_create_options$allows_duplicates(),
        offsets_filter_list = tiledb::tiledb_filter_list(tiledb_create_options$offsets_filters()),
        validity_filter_list = tiledb::tiledb_filter_list(tiledb_create_options$validity_filters())
        ## enumerations = if (any(!sapply(levels, is.null))) levels else NULL
        )

      for (field_name in attr_column_names) {
          fieldtype <- schema$GetFieldByName(field_name)$type
          if (is(fieldtype, "DictionaryType")) {
              tiledb::tiledb_array_schema_set_enumeration_empty(schema = tdb_schema,
                                                                attr = tdb_attrs[[field_name]],
                                                                enum_name = field_name,
                                                                type_str = "UTF8",
                                                                ordered = fieldtype$ordered)
          }
      }

      # create array
      tiledb::tiledb_array_create(uri = self$uri, schema = tdb_schema)
      self$open("WRITE", internal_use_only = "allowed_use")
      private$write_object_type_metadata()
      self
    },

    #' @description Write (lifecycle: experimental)
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

      df <- as.data.frame(values)[schema_names]
      arr <- self$object

      has_enums <- tiledb::tiledb_array_has_enumeration(arr)
      if (any(has_enums)) {       # if enumerations exists in array
          attrs <- tiledb::attrs(tiledb::schema(arr))
          if (!tiledb::tiledb_array_is_open(arr)) arr <- tiledb::tiledb_array_open(arr, "READ")
          for (attr_name in names(attrs)) {
              if (has_enums[attr_name]) {
                  old_enum <- tiledb::tiledb_attribute_get_enumeration(attrs[[attr_name]], arr)
                  new_enum <- levels(values[[attr_name]]$as_vector())
                  added_enum <- setdiff(new_enum, old_enum)
                  if (length(added_enum) > 0) {
                      ase <- tiledb::tiledb_array_schema_evolution()
                      ase <- tiledb::tiledb_array_schema_evolution_extend_enumeration(ase, arr, attr_name, added_enum)
                      tiledb::tiledb_array_schema_evolution_array_evolve(ase, self$uri)
                      df[, attr_name] <- factor(df[, attr_name], levels = unique(c(old_enum,new_enum)))
                  }
              }
          }
      }

      arr[] <- df
      # tiledb-r always closes the array after a write operation so we need to
      # manually reopen it until close-on-write is optional
      self$open("WRITE", internal_use_only = "allowed_use")
      invisible(self)
    },

    #' @description Read (lifecycle: experimental)
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

      cfg <- as.character(tiledb::config(self$tiledbsoma_ctx$context()))
      rl <- sr_setup(uri = self$uri,
                     config = cfg,
                     colnames = column_names,
                     qc = value_filter,
                     dim_points = coords,
                     timestamp_end = private$tiledb_timestamp,
                     loglevel = log_level)
      private$ctx_ptr <- rl$ctx
      TableReadIter$new(rl$sr)
    },

    #' @description Update (lifecycle: experimental)
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
      private$reopen(mode = "READ")
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
      private$reopen(mode = "WRITE")
      spdl::info("[SOMADataFrame update]: Writing new data")
      self$write(values)
    }

  ),

  private = list(

    # @description Validate schema (lifecycle: experimental)
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
