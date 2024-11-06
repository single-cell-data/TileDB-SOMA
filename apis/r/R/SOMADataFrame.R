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
    #' @param domain An optional list of 2-element vectors specifying the domain of each index
    #' column. Each vector should be a pair consisting of the minimum and maximum values storable in
    #' the index column. For example, if there is a single int64-valued index column, then `domain`
    #' might be `c(100, 200)` to indicate that values between 100 and 200, inclusive, can be stored
    #' in that column.  If provided, this list must have the same length as `index_column_names`,
    #' and the index-column domain will be as specified.  If omitted entirely, or if `NULL` in a given
    #' dimension, the corresponding index-column domain will use the minimum and maximum possible
    #' values for the column's datatype.  This makes a `DataFrame` growable.
    #' @template param-platform-config
    #' @param internal_use_only Character value to signal this is a 'permitted' call,
    #' as `create()` is considered internal and should not be called directly.
    create = function(
      schema,
      index_column_names = c("soma_joinid"),
      domain = NULL,
      platform_config = NULL,
      internal_use_only = NULL
    ) {
      if (is.null(internal_use_only) || internal_use_only != "allowed_use") {
        stop(paste(
          "Use of the create() method is for internal use only. Consider using a",
          "factory method as e.g. 'SOMADataFrameCreate()'."
        ), call. = FALSE)
      }
      schema <- private$validate_schema(schema, index_column_names)

      stopifnot(
        "domain must be NULL or a named list, with values being 2-element vectors or NULL" = is.null(domain) ||
          ( # Check that `domain` is a list of length `length(index_column_names)`
            # where all values are named after `index_column_names`
            # and all values are `NULL` or a two-length atomic non-factor vector
            rlang::is_list(domain, n = length(index_column_names)) &&
              identical(sort(names(domain)), sort(index_column_names)) &&
              all(vapply_lgl(
                domain,
                function(x) is.null(x) || (is.atomic(x) && !is.factor(x) && length(x) == 2L)
              ))
          )
      )

      attr_column_names <- setdiff(schema$names, index_column_names)
      stopifnot(
        "At least one non-index column must be defined in the schema" =
          length(attr_column_names) > 0
      )

      if ("soma_joinid" %in% index_column_names && !is.null(domain)) {
        lower <- domain[["soma_joinid"]][1]
        stopifnot("The lower bound for soma_joinid domain must be 0" = lower == 0)
      }

      # Parse the tiledb/create/ subkeys of the platform_config into a handy,
      # typed, queryable data structure.
      tiledb_create_options <- TileDBCreateOptions$new(platform_config)

      # We currently pass domain and extent values in an arrow table (i.e. data.frame alike)
      # where each dimension is one column (of the same type as in the schema followed by:
      # * Before the new shape feature: three values for the domain pair and the extent;
      # * After the new shape feature: five values for the maxdomain pair, extent, and domain.
      dom_ext_tbl <- get_domain_and_extent_dataframe(
        schema,
        ind_col_names = index_column_names,
        domain = domain,
        tdco = tiledb_create_options
      )

      ## we transfer to the arrow table via a pair of array and schema pointers
      dnaap <- nanoarrow::nanoarrow_allocate_array()
      dnasp <- nanoarrow::nanoarrow_allocate_schema()
      arrow::as_record_batch(dom_ext_tbl)$export_to_c(dnaap, dnasp)

      ## we need a schema pointer to transfer the schema information
      nasp <- nanoarrow::nanoarrow_allocate_schema()
      schema$export_to_c(nasp)

      spdl::debug("[SOMADataFrame$create] about to create schema from arrow")
      createSchemaFromArrow(
        uri = self$uri,
        nasp = nasp,
        nadimap = dnaap,
        nadimsp = dnasp,
        sparse = TRUE,
        datatype = "SOMADataFrame",
        pclst = tiledb_create_options$to_list(FALSE),
        ctxxp = private$.soma_context,
        tsvec = self$.tiledb_timestamp_range
      )

      spdl::debug("[SOMADataFrame$create] about to call write_object_type_metadata")
      private$write_object_type_metadata()

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
        ctxxp = private$.soma_context,
        arraytype = "SOMADataFrame",
        config = NULL,
        tsvec = self$.tiledb_timestamp_range
      )

      invisible(self)
    },

    #' @description Read (lifecycle: maturing)
    #' Read a user-defined subset of data, addressed by the dataframe indexing
    #' column, and optionally filtered.
    #' @param coords Optional named list of indices specifying the rows to read;
    #' each (named) list element corresponds to a dimension of the same name.
    #' @param column_names Optional character vector of column names to return.
    #' @param value_filter Optional string containing a logical expression that is used
    #' to filter the returned values. See [`tiledb::parse_query_condition`] for
    #' more information.
    #' @template param-result-order
    #' @param iterated Option boolean indicated whether data is read in call (when
    #' `FALSE`, the default value) or in several iterated steps.
    #' @param log_level Optional logging level with default value of `"warn"`.
    #' @return arrow::\link[arrow]{Table} or \link{TableReadIter}
    read = function(
      coords = NULL,
      column_names = NULL,
      value_filter = NULL,
      result_order = "auto",
      iterated = FALSE,
      log_level = "auto"
    ) {
      private$check_open_for_read()

      result_order <- match_query_layout(result_order)

      ## if unnamed set names
      if (!is.null(coords)) {
        if (!is.list(coords)) {
          coords <- list(coords)
        }
        if (is.null(names(coords))) {
          names(coords) <- self$dimnames()
        }
      }

      stopifnot(
        ## check columns
        "'column_names' must only contain valid dimension or attribute columns" =
          is.null(column_names) || all(column_names %in% c(self$dimnames(), self$attrnames()))
      )

      coords <- validate_read_coords(coords, dimnames = self$dimnames(), schema = self$schema())

      if (!is.null(value_filter)) {
        value_filter <- validate_read_value_filter(value_filter)
        parsed <- do.call(
          what = parse_query_condition,
          args = list(expr = value_filter, schema = self$schema(), somactx = private$.soma_context)
        )
        value_filter <- parsed@ptr
      }
      spdl::debug(
        "[SOMADataFrame$read] calling sr_setup for {} at ({},{})", self$uri,
        private$tiledb_timestamp[1], private$tiledb_timestamp[2]
      )
      sr <- sr_setup(
        uri = self$uri,
        private$.soma_context,
        colnames = column_names,
        qc = value_filter,
        dim_points = coords,
        timestamprange = self$.tiledb_timestamp_range, # NULL or two-elem vector
        loglevel = log_level
      )
      TableReadIter$new(sr)
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

      # Leave state unmodified
      # TODO: this issue will automatically go away on https://github.com/single-cell-data/TileDB-SOMA/issues/3059
      omode <- self$mode()
      on.exit(self$reopen(mode = omode))

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

      drop_cols_for_clib <- drop_cols
      add_cols_types_for_clib <-
        add_cols_enum_value_types_for_clib <-
        add_cols_enum_ordered_for_clib <- vector("list", length = length(add_cols))
      names(add_cols_types_for_clib) <-
        names(add_cols_enum_value_types_for_clib) <-
        names(add_cols_enum_ordered_for_clib) <- add_cols

      # Add columns
      for (add_col in add_cols) {
        col_type <- new_schema$GetFieldByName(add_col)$type

        if (inherits(col_type, "DictionaryType")) {
          spdl::debug(
            "[SOMADataFrame update]: adding enum column '{}' index type '{}' value type '{}' ordered {}",
            add_col, col_type$index_type$name, col_type$value_type$name, col_type$ordered
          )

          add_cols_types_for_clib[[add_col]] <- col_type$index_type$name
          add_cols_enum_value_types_for_clib[[add_col]] <- col_type$value_type$name
          add_cols_enum_ordered_for_clib[[add_col]] <- col_type$ordered
        } else {
          spdl::debug("[SOMADataFrame update]: adding column '{}' type '{}'", add_col, col_type$name)

          add_cols_types_for_clib[[add_col]] <- col_type$name
        }
      }

      if (length(drop_cols_for_clib) > 0 || length(add_cols_types_for_clib) > 0) {
        c_update_dataframe_schema(
          self$uri,
          private$.soma_context,
          drop_cols_for_clib,
          Filter(Negate(is.null), add_cols_types_for_clib),
          Filter(Negate(is.null), add_cols_enum_value_types_for_clib),
          Filter(Negate(is.null), add_cols_enum_ordered_for_clib)
        )
      }

      # Reopen array for writing with new schema
      self$reopen(mode = "WRITE")

      spdl::debug("[SOMADataFrame update]: Writing new data")
      self$write(values)
    },

    #' @description Retrieve the shape; as \code{SOMADataFrames} are shapeless,
    #' simply raises an error
    #'
    #' @return None, instead a \code{\link{.NotYetImplemented}()} error is raised
    #'
    shape = function() {
      stop(errorCondition(
        "'SOMADataFrame$shape()' is not implemented yet",
        class = "notYetImplementedError"
      ))
    },

    #' @description Retrieve the maxshape; as \code{SOMADataFrames} are shapeless,
    #' simply raises an error
    #'
    #' @return None, instead a \code{\link{.NotYetImplemented}()} error is raised
    #'
    maxshape = function() {
      stop(errorCondition(
        "'SOMADataFrame$maxshape()' is not implemented",
        class = "notYetImplementedError"
      ))
    },

    #' @description Returns a named list of minimum/maximum pairs, one per index
    #' column, currently storable on each index column of the dataframe. These
    #' can be resized up to `maxdomain`.
    #' (lifecycle: maturing)
    #' @return Named list of minimum/maximum values.
    domain = function() {
      as.list(
        arrow::as_record_batch(
          arrow::as_arrow_table(
            domain(self$uri, private$.soma_context)
          )
        )
      )
    },

    #' @description Returns a named list of minimum/maximum pairs, one per index
    #' column, which are the limits up to which the dataframe can have its
    #' domain resized.
    #' (lifecycle: maturing)
    #' @return Named list of minimum/maximum values.
    maxdomain = function() {
      as.list(
        arrow::as_record_batch(
          arrow::as_arrow_table(
            maxdomain(self$uri, private$.soma_context)
          )
        )
      )
    },

    #' @description Returns TRUE if the array has the upgraded resizeable domain
    #' feature from TileDB-SOMA 1.15: the array was created with this support,
    #' or it has had ``upgrade_domain`` applied to it.
    #' (lifecycle: maturing)
    #' @return Logical
    tiledbsoma_has_upgraded_domain = function() {
      has_current_domain(self$uri, private$.soma_context)
    },

    #' @description Increases the shape of the dataframe on the ``soma_joinid``
    #' index column, if it indeed is an index column, leaving all other index
    #' columns as-is. If the ``soma_joinid`` is not an index column, no change is
    #' made.  This is a special case of ``upgrade_domain`` (WIP for 1.15), but
    #' simpler to keystroke, and handles the most common case for dataframe
    #' domain expansion.  Raises an error if the dataframe doesn't already have a
    #' domain: in that case please call ``tiledbsoma_upgrade_domain`` (WIP for
    #' 1.15).
    #' @param new_shape An integer, greater than or equal to 1 + the
    #' `soma_joinid` domain slot.
    #' @return No return value
    tiledbsoma_resize_soma_joinid_shape = function(new_shape) {
      stopifnot("'new_shape' must be an integer" = rlang::is_integerish(new_shape, n = 1) ||
        (bit64::is.integer64(new_shape) && length(new_shape) == 1))
      # Checking slotwise new shape >= old shape, and <= max_shape, is already done in libtiledbsoma
      invisible(
        resize_soma_joinid_shape(
          self$uri, new_shape, .name_of_function(), private$.soma_context))
    },

    #' @description Allows you to set the domain of a `SOMADataFrame`, when the
    #' `SOMADataFrame` does not have a domain set yet.  The argument must be a
    #' tuple of pairs of low/high values for the desired domain, one pair per
    #' index column. For string index columns, you must offer the low/high pair
    #' as `("", "")`, or as `NULL`.  If ``check_only`` is ``True``, returns
    #' whether the operation would succeed if attempted, and a reason why it
    #' would not. The domain being requested must be contained within what
    #' `maxdomain` returns.
    #' @param new_domain A named list, keyed by index-column name, with values
    #' being two-element vectors containing the desired lower and upper bounds
    #' for the domain.
    #' @return No return value
    tiledbsoma_upgrade_domain = function(new_domain) {
      # stopifnot("'new_domain' must be CODE ME UP PLZ" = ...
      # Checking slotwise new shape >= old shape, and <= max_shape, is already
      # done in libtiledbsoma

      pyarrow_domain_table <- private$upgrade_or_change_domain_helper(
        new_domain, "tiledbsoma_upgrade_domain"
      )

      invisible(
        upgrade_or_change_domain(
          self$uri,
          FALSE,
          pyarrow_domain_table$array,
          pyarrow_domain_table$schema,
          .name_of_function(),
          private$.soma_context
        )
      )
    },

    #' @description Allows you to set the domain of a `SOMADataFrame`, when the
    #' `SOMADataFrame` already has a domain set yet.  The argument must be a
    #' tuple of pairs of low/high values for the desired domain, one pair per
    #' index column. For string index columns, you must offer the low/high pair
    #' as `("", "")`, or as `NULL`.  If ``check_only`` is ``True``, returns
    #' whether the operation would succeed if attempted, and a reason why it
    #' would not. The return value from `domain` must be contained within
    #' the requested `new_domain`, and the requested `new_domain` must be
    #' contained within the return value from `maxdomain`.
    #' @param new_domain A named list, keyed by index-column name, with values
    #' being two-element vectors containing the desired lower and upper bounds
    #' for the domain.
    #' @return No return value
    change_domain = function(new_domain) {
      # stopifnot("'new_domain' must be CODE ME UP PLZ" = ...
      # Checking slotwise new shape >= old shape, and <= max_shape, is already
      # done in libtiledbsoma

      pyarrow_domain_table <- private$upgrade_or_change_domain_helper(
        new_domain, tiledbsoma_upgrade_domain
      )

      invisible(
        upgrade_or_change_domain(
          self$uri,
          TRUE,
          pyarrow_domain_table$array,
          pyarrow_domain_table$schema,
          .name_of_function(),
          private$.soma_context
        )
      )
    }
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
          all(!startsWith(setdiff(schema$names, "soma_joinid"), "soma_")) ||
          isTRUE(getOption("tiledbsoma.write_soma.internal", default = FALSE))
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

    # Converts the user-level tuple of low/high pairs into a pyarrow table
    # suitable for calling libtiledbsoma.
    upgrade_or_change_domain_helper = function(
      new_domain,
      function_name_for_messages
    ) {
      dimnames <- self$dimnames()

      # Check user-provided domain against dataframe domain.
      stopifnot(
        "new_domain must be a named list, with values being 2-element vectors or NULL, with names the same as the dataframe's index-column names" =
        rlang::is_list(new_domain, n = length(dimnames)) &&
          identical(sort(names(new_domain)), sort(dimnames)) &&
          all(vapply_lgl(
            new_domain,
            function(x) is.null(x) || (is.atomic(x) && !is.factor(x) && length(x) == 2L)
          ))
      )

      # From the dataframe's schema, extract the subschema for only index
      # columns (TileDB dimensions).
      full_schema <- self$schema()
      dim_schema_list <- list()
      for (dimname in dimnames) {
        dim_schema_list[[dimname]] <- full_schema[[dimname]]
      }
      dim_schema <- arrow::schema(dim_schema_list)

      # Get the user's new_domain list with keys in the same order as dim_names.
      ordered_new_domain = list()
      for (dimname in dimnames) {
        # * Domain cannot be specified for string-type index columns.
        # * So we let them say `NULL` rather than `c("", "")`.
        # * But R list semantics are `mylist[[key]] <- NULL` results in nothing
        #   being set at that key.
        if (is.null(new_domain[[dimname]]) && full_schema[[dimname]]$type$ToString() %in% c("string", "large_string", "utf8", "large_utf8")) {
          ordered_new_domain[[dimname]] <- c("", "")
        } else {
          ordered_new_domain[[dimname]] <- new_domain[[dimname]]
        }
      }

      pyarrow_table <- arrow::arrow_table(as.data.frame(ordered_new_domain), schema=dim_schema)

      # We transfer to the arrow table via a pair of array and schema pointers
      dnaap <- nanoarrow::nanoarrow_allocate_array()
      dnasp <- nanoarrow::nanoarrow_allocate_schema()
      arrow::as_record_batch(pyarrow_table)$export_to_c(dnaap, dnasp)

      return(list(array=dnaap, schema=dnasp))
    }
  )
)
