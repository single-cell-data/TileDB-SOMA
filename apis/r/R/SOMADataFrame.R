#' SOMADataFrame
#'
#' @description A SOMA data frame  is a multi-column table that must contain a
#' column called \dQuote{\code{soma_joinid}} of type \code{int64}, which
#' contains a unique value for each row and is intended to act as a join key for
#' other objects, such as \code{\link{SOMASparseNDArray}} (lifecycle: maturing).
#'
#' @param new_domain A named list, keyed by index-column name, with values
#' being two-element vectors containing the desired lower and upper bounds
#' for the domain.
#' @param check_only If true, does not apply the operation, but only reports
#' whether it would have succeeded.
#'
#' @export
#'
#' @inherit SOMADataFrameCreate examples
#'
SOMADataFrame <- R6::R6Class(
  classname = "SOMADataFrame",
  inherit = SOMAArrayBase,
  public = list(
    #' @description Create a SOMA data frame (lifecycle: maturing).\cr
    #' \cr
    #' \strong{Note}: \code{$create()} is considered internal and should not be
    #' called directly; use factory functions
    #' (eg. \code{\link{SOMADataFrameCreate}()}) instead.
    #'
    #' @param schema An \link[arrow:schema]{Arrow schema}.
    #' @param index_column_names A vector of column names to use as user-defined
    #' index columns. All named columns must exist in the schema, and at least
    #' one index column name is required.
    #' @param domain An optional list specifying the domain of each
    #' index column. Each slot in the list must have its name being the name
    #' of an index column, and its value being be a length-two vector
    #' consisting of the minimum and maximum values storable in the index
    #' column. For example, if there is a single int64-valued index column
    #' \code{soma_joinid}, then \code{domain} might be
    #' \code{list(soma_joinid=c(100, 200))} to indicate that values between 100
    #' and 200, inclusive, can be stored in that column.  If provided, this
    #' sequence must have the same length as \code{index_column_names}, and the
    #' index-column domain will be as specified. Omitting or setting the domain
    #' to \code{NULL} is deprecated. See also \code{change_domain} which allows
    #' you to expand the domain after create.
    #' @template param-platform-config
    #'
    #' @return Returns \code{self}.
    #'
    create = function(
      schema,
      index_column_names = c("soma_joinid"),
      domain = NULL,
      platform_config = NULL
    ) {
      envs <- unique(vapply(
        X = unique(sys.parents()),
        FUN = function(n) environmentName(environment(sys.function(n))),
        FUN.VALUE = character(1L)
      ))
      if (!"tiledbsoma" %in% envs) {
        stop(
          paste(
            strwrap(private$.internal_use_only("create", "collection")),
            collapse = '\n'
          ),
          call. = FALSE
        )
      }

      schema <- private$validate_schema(schema, index_column_names)

      if (is.null(domain)) {
        .deprecate(
          when = "2.1.0",
          what = "create(domain = 'must be a named list')",
        )
      }
      if (!(is.null(domain) || .is_domain(domain, index_column_names))) {
        stop(
          "domain must be NULL or a named list, with values being 2-element vectors or NULL"
        )
      }

      attr_column_names <- setdiff(schema$names, index_column_names)
      stopifnot(
        "At least one non-index column must be defined in the schema" = length(
          attr_column_names
        ) >
          0
      )

      if ("soma_joinid" %in% index_column_names && !is.null(domain)) {
        lower_bound <- domain[["soma_joinid"]][1]
        upper_bound <- domain[["soma_joinid"]][2]
        stopifnot(
          "The lower bound for soma_joinid domain must be >= 0" = lower_bound >=
            0,
          "The upper bound for soma_joinid domain must be >= 0" = upper_bound >=
            0,
          "The upper bound for soma_joinid domain must be >= the lower bound" = upper_bound >=
            lower_bound
        )
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

      soma_debug("[SOMADataFrame$create] about to create schema from arrow")
      createSchemaFromArrow(
        uri = self$uri,
        nasp = nasp,
        nadimap = dnaap,
        nadimsp = dnasp,
        sparse = TRUE,
        datatype = "SOMADataFrame",
        pclst = tiledb_create_options$to_list(FALSE),
        ctxxp = private$.context$handle,
        tsvec = self$.tiledb_timestamp_range
      )

      soma_debug(
        "[SOMADataFrame$create] about to call write_object_type_metadata"
      )
      self$open("WRITE")
      private$write_object_type_metadata()
      self$reopen("WRITE", tiledb_timestamp = self$tiledb_timestamp)

      return(self)
    },

    #' @description Write values to the data frame (lifecycle: maturing).
    #'
    #' @param values An \link[arrow:Table]{Arrow table} or
    #' \link[arrow:RecordBatch]{Arrow record batch} containing all columns,
    #' including any index columns. The schema for \code{values} must match the
    #' schema for the data frame.
    #'
    #' @return Invisibly returns \code{self}.
    #'
    write = function(values) {
      private$.check_handle()

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
        "'values' must be an Arrow Table or RecordBatch" = (is_arrow_table(
          values
        ) ||
          is_arrow_record_batch(values)),
        "All columns in 'values' must be defined in the schema" = all(
          col_names %in% schema_names
        ),
        "All schema fields must be present in 'values'" = all(
          schema_names %in% col_names
        )
      )

      ## we transfer to the arrow table via a pair of array and schema pointers
      naap <- nanoarrow::nanoarrow_allocate_array()
      nasp <- nanoarrow::nanoarrow_allocate_schema()
      arrow::as_record_batch(values)$export_to_c(naap, nasp)
      writeArrayFromArrow(
        soma_array = private$.handle,
        naap = naap,
        nasp = nasp,
        arraytype = "SOMADataFrame"
      )

      return(invisible(self))
    },

    #' @description Read a user-defined subset of data, addressed by the
    #' data frame indexing column, and optionally filtered
    #' (lifecycle: maturing).
    #'
    #' @param coords Optional named list of indices specifying the rows to read;
    #' each (named) list element corresponds to a dimension of the same name.
    #' @param column_names Optional character vector of column names to return.
    #' @param value_filter Optional string containing a logical expression that
    #' is used to filter the returned values. See
    #' \code{\link[tiledb:parse_query_condition]{tiledb::parse_query_condition}()}
    #' for more information.
    #' @template param-result-order
    #' @param log_level Optional logging level with default value of
    #' \dQuote{\code{warn}}.
    #'
    #' @return An \link[arrow:Table]{Arrow table} or \code{\link{TableReadIter}}
    #'
    read = function(
      coords = NULL,
      column_names = NULL,
      value_filter = NULL,
      result_order = "auto",
      log_level = "auto"
    ) {
      private$.check_open_for_read()

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
        "'column_names' must only contain valid dimension or attribute columns" = is.null(
          column_names
        ) ||
          all(column_names %in% c(self$dimnames(), self$attrnames()))
      )

      coords <- validate_read_coords(
        coords,
        dimnames = self$dimnames(),
        schema = self$schema()
      )

      if (!is.null(value_filter)) {
        value_filter <- validate_read_value_filter(value_filter)
        parsed <- do.call(
          what = parse_query_condition,
          args = list(
            expr = value_filter,
            schema = self$schema(),
            somactx = private$.context$handle
          )
        )
        value_filter <- parsed@ptr
      }

      if (is.null(self$.tiledb_timestamp_range)) {
        soma_debug(sprintf(
          "[SOMADataFrame$read] calling mq_setup for %s",
          self$uri
        ))
      } else {
        soma_debug(sprintf(
          "[SOMADataFrame$read] calling mq_setup for %s at (%s, %s)",
          self$uri,
          self$.tiledb_timestamp_range[1],
          self$.tiledb_timestamp_range[2]
        ))
      }

      sr <- mq_setup(
        uri = self$uri,
        private$.context$handle,
        colnames = column_names,
        qc = value_filter,
        dim_points = coords,
        timestamprange = self$.tiledb_timestamp_range, # NULL or two-elem vector
        loglevel = log_level
      )
      return(TableReadIter$new(sr))
    },

    #' @description Update (lifecycle: maturing).
    #'
    #' @details
    #' Update the existing \code{SOMADataFrame} to add or remove columns based
    #' on the input:
    #' \itemize{
    #'  \item columns present in the current the \code{SOMADataFrame} but absent
    #'   from the new \code{values} will be dropped.
    #'  \item columns absent in current \code{SOMADataFrame} but present in the
    #'   new \code{values} will be added.
    #'  \item any columns present in both will be left alone, with the
    #'   exception that if \code{values} has a different type for the column,
    #'   the entire update will fail because attribute types cannot be changed.
    #' }
    #' Furthermore, \code{values} must contain the same number of rows as the
    #' current \code{SOMADataFrame}.
    #'
    #' @param values A data frame, \link[arrow:Table]{Arrow table}, or
    #' \link[arrow:RecordBatch]{Arrow record batch}.
    #' @param row_index_name An optional scalar character. If provided, and if
    #' the \code{values} argument is a data frame with row names, then the row
    #' names will be extracted and added as a new column to the data frame
    #' prior to performing the update. The name of this new column will be set
    #' to the value specified by \code{row_index_name}.
    #'
    #' @return Invisibly returns \code{NULL}
    #'
    update = function(values, row_index_name = NULL) {
      private$.check_open_for_write()
      stopifnot(
        "'values' must be a data.frame, Arrow Table or RecordBatch" = is.data.frame(
          values
        ) ||
          is_arrow_table(values) ||
          is_arrow_record_batch(values)
      )

      # Leave state unmodified
      # TODO: this issue will automatically go away on https://github.com/single-cell-data/TileDB-SOMA/issues/3059
      omode <- self$mode()
      on.exit(self$reopen(mode = omode))

      if (is.data.frame(values)) {
        if (!is.null(row_index_name)) {
          stopifnot(
            "'row_index_name' must be NULL or a scalar character vector" = is_scalar_character(
              row_index_name
            ),
            "'row_index_name' conflicts with an existing column name" = !row_index_name %in%
              colnames(values)
          )
          values[[row_index_name]] <- rownames(values)
        }
        values <- arrow::as_arrow_table(values)
      }

      # Retrieve existing soma_joinids from array to:
      # - validate number of rows in values matches number of rows in array
      # - add original soma_joinids to values if not present
      soma_debug("[SOMADataFrame update]: Retrieving existing soma_joinids")
      self$reopen(mode = "READ")
      joinids <- self$read(column_names = "soma_joinid")$concat()$soma_joinid
      if (length(joinids) != nrow(values)) {
        stop(
          "Number of rows in 'values' must match number of rows in array:\n",
          "  - Number of rows in array: ",
          length(joinids),
          "\n",
          "  - Number of rows in 'values': ",
          nrow(values),
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
      add_cols_types_for_clib <- stats::setNames(
        vector("list", length(add_cols)),
        nm = add_cols
      )
      add_cols_enum_value_types_for_clib <- stats::setNames(
        vector("list", length(add_cols)),
        nm = add_cols
      )
      add_cols_enum_ordered_for_clib <- stats::setNames(
        vector("list", length(add_cols)),
        nm = add_cols
      )

      # Add columns
      for (add_col in add_cols) {
        col_type <- new_schema$GetFieldByName(add_col)$type

        if (inherits(col_type, "DictionaryType")) {
          soma_debug(sprintf(
            "[SOMADataFrame update]: adding enum column '%s' index type '%s' value type '%s' ordered %s",
            add_col,
            col_type$index_type$name,
            col_type$value_type$name,
            col_type$ordered
          ))

          add_cols_types_for_clib[[add_col]] <- col_type$index_type$name
          add_cols_enum_value_types_for_clib[[
            add_col
          ]] <- col_type$value_type$name
          add_cols_enum_ordered_for_clib[[add_col]] <- col_type$ordered
        } else {
          soma_debug(sprintf(
            "[SOMADataFrame update]: adding column '%s' type '%s'",
            add_col,
            col_type$name
          ))

          add_cols_types_for_clib[[add_col]] <- col_type$name
        }
      }

      if (
        length(drop_cols_for_clib) > 0 || length(add_cols_types_for_clib) > 0
      ) {
        c_update_dataframe_schema(
          self$uri,
          private$.context$handle,
          drop_cols_for_clib,
          Filter(Negate(is.null), add_cols_types_for_clib),
          Filter(Negate(is.null), add_cols_enum_value_types_for_clib),
          Filter(Negate(is.null), add_cols_enum_ordered_for_clib)
        )
      }

      # Reopen array for writing with new schema
      self$reopen(mode = "WRITE")

      soma_debug("[SOMADataFrame update]: Writing new data")
      self$write(values)
    },

    #' @description Get the levels for an enumerated (\code{factor}) column.
    #'
    #' @param column_names Optional character vector of column names to pull
    #' enumeration levels for; defaults to all enumerated columns.
    #' @param simplify Simplify the result down to a vector or matrix.
    #'
    #' @return If \code{simplify} returns one of the following:
    #' \itemize{
    #'  \item a vector of there is only one enumerated column.
    #'  \item a matrix if there are multiple enumerated columns with the same
    #'   number of levels.
    #'  \item a named list if there are multiple enumerated columns with
    #'   differing numbers of levels.
    #' }
    #' Otherwise, returns a named list.
    #'
    levels = function(column_names = NULL, simplify = TRUE) {
      stopifnot(
        "'simplify' must be TRUE or FALSE" = isTRUE(simplify) ||
          isFALSE(simplify)
      )
      enums <- c_attributes_enumerated(self$uri, private$.context$handle)
      if (!any(enums)) {
        rlang::warn("No enumerated columns present")
        return(NULL)
      }
      factors <- names(enums[enums])
      column_names <- column_names %||% factors
      column_names <- rlang::arg_match(
        column_names,
        values = factors,
        multiple = TRUE
      )
      levels <- sapply(
        X = column_names,
        FUN = c_attribute_enumeration_levels,
        uri = self$uri,
        ctxxp = private$.context$handle,
        simplify = FALSE,
        USE.NAMES = TRUE
      )
      if (simplify) {
        if (length(levels) == 1L) {
          return(levels[[1L]])
        }
        return(simplify2array(levels, higher = FALSE))
      }
      return(levels)
    },

    #' @description Retrieve the shape; as \code{SOMADataFrames} are shapeless,
    #' simply raises an error.
    #'
    #' @return None, instead a \code{\link{.NotYetImplemented}()} error
    #' is raised.
    #'
    shape = function() {
      stop(errorCondition(
        "'SOMADataFrame$shape()' is not implemented yet",
        class = "notYetImplementedError"
      ))
    },

    #' @description Retrieve the max shape; as \code{SOMADataFrames} are
    #' shapeless, simply raises an error.
    #'
    #' @return None, instead a \code{\link{.NotYetImplemented}()} error
    #' is raised.
    #'
    maxshape = function() {
      stop(errorCondition(
        "'SOMADataFrame$maxshape()' is not implemented",
        class = "notYetImplementedError"
      ))
    },

    #' @description Returns a named list of minimum/maximum pairs, one per index
    #' column, currently storable on each index column of the data frame. These
    #' can be resized up to \code{maxdomain} (lifecycle: maturing).
    #'
    #' @return Named list of minimum/maximum values.
    #'
    domain = function() {
      private$.check_handle()
      return(as.list(arrow::as_record_batch(arrow::as_arrow_table(domain(
        private$.handle
      )))))
    },

    #' @description Returns a named list of minimum/maximum pairs, one per index
    #' column, which are the limits up to which the data frame can have its
    #' domain resized (lifecycle: maturing).
    #'
    #' @return Named list of minimum/maximum values.
    #'
    maxdomain = function() {
      private$.check_handle()
      return(as.list(arrow::as_record_batch(arrow::as_arrow_table(maxdomain(
        private$.handle
      )))))
    },

    #' @description Test if the array has the upgraded resizeable domain feature
    #' from TileDB-SOMA 1.15, the array was created with this support, or it has
    #' had \code{$upgrade_domain()} applied to it (lifecycle: maturing).
    #'
    #' @return Returns \code{TRUE} if the array has the upgraded resizable
    #' domain feature; otherwise, returns \code{FALSE}.
    #'
    tiledbsoma_has_upgraded_domain = function() {
      private$.check_handle()
      has_current_domain(private$.handle)
    },

    #' @description Increases the shape of the data frame on the
    #' \code{soma_joinid} index column, if it indeed is an index column, leaving
    #' all other index columns as-is. If the \code{soma_joinid} is not an index
    #' column, no change is made. This is a special case of
    #' \code{upgrade_domain()}, but simpler to keystroke, and
    #' handles the most common case for data frame domain expansion. Raises an
    #' error if the data frame doesn't already have a domain; in that case
    #' please call \code{$tiledbsoma_upgrade_domain()}.
    #'
    #' @param new_shape An integer, greater than or equal to 1 + the
    #' \code{soma_joinid} domain slot.
    #'
    #' @return Invisibly returns \code{NULL}
    #'
    tiledbsoma_resize_soma_joinid_shape = function(new_shape) {
      private$.check_handle()
      stopifnot(
        "'new_shape' must be an integer" = rlang::is_integerish(
          new_shape,
          n = 1
        ) ||
          (bit64::is.integer64(new_shape) && length(new_shape) == 1)
      )
      # Checking slotwise new shape >= old shape, and <= max_shape, is already done in libtiledbsoma
      return(invisible(resize_soma_joinid_shape(
        private$.handle,
        new_shape,
        .name_of_function()
      )))
    },

    #' @description Allows you to set the domain of a \code{SOMADataFrame},
    #' when the \code{SOMADataFrame} does not have a domain set yet. The
    #' argument must be a list of pairs of low/high values for the desired
    #' domain, one pair per index column. For string index columns, you must
    #' offer the low/high pair as \code{c("", "")}, or as \code{NULL}. If
    #' \code{check_only} is \code{True}, returns whether the operation would
    #' succeed if attempted, or a reason why it would not. The domain being
    #' requested must be contained within what \code{$maxdomain()} returns.
    #'
    #' @return If \code{check_only}, returns the empty string if no error is
    #' detected, else a description of the error. Otherwise, invisibly returns
    #' \code{NULL}
    #'
    tiledbsoma_upgrade_domain = function(new_domain, check_only = FALSE) {
      private$.check_handle()

      pyarrow_domain_table <- private$upgrade_or_change_domain_helper(
        new_domain,
        "tiledbsoma_upgrade_domain"
      )

      reason_string <- upgrade_or_change_domain(
        private$.handle,
        FALSE,
        pyarrow_domain_table$array,
        pyarrow_domain_table$schema,
        .name_of_function(),
        check_only
      )

      if (isTRUE(check_only)) {
        return(reason_string)
      }

      # The return value from upgrade_or_change_domain without check_only is
      # always "", or it raises an error trying.
      return(invisible(NULL))
    },

    #' @description Allows you to set the domain of a \code{SOMADataFrame}, when
    #' the \code{SOMADataFrame} already has a domain set yet. The argument must
    #' be a list of pairs of low/high values for the desired domain, one pair
    #' per index column. For string index columns, you must offer the low/high
    #' pair as \code{c("", "")}, or as \code{NULL}. If \code{check_only} is
    #' \code{True}, returns whether the operation would succeed if attempted,
    #' or a reason why it would not. The return value from \code{domain} must be
    #' contained within the requested \code{new_domain}, and the requested
    #' \code{new_domain} must be contained within the return value from
    #' \code{$maxdomain()} (lifecycle: maturing).
    #'
    #' @return If \code{check_only}, returns the empty string if no error is
    #' detected, else a description of the error. Otherwise, invisibly returns
    #' \code{NULL}
    #'
    change_domain = function(new_domain, check_only = FALSE) {
      private$.check_handle()

      pyarrow_domain_table <- private$upgrade_or_change_domain_helper(
        new_domain,
        tiledbsoma_upgrade_domain
      )

      reason_string <- upgrade_or_change_domain(
        private$.handle,
        TRUE,
        pyarrow_domain_table$array,
        pyarrow_domain_table$schema,
        .name_of_function(),
        check_only
      )

      if (isTRUE(check_only)) {
        return(reason_string)
      }

      # The return value from upgrade_or_change_domain without check_only is
      # always "", or it raises an error trying.
      return(invisible(NULL))
    }
  ),
  private = list(
    # @description Open the handle for the C++ interface
    .open_handle = function(open_mode, timestamp) {
      private$.handle <- open_dataframe_handle(
        self$uri,
        open_mode,
        private$.context$handle,
        timestamp
      )
    },

    # @description Validate schema (lifecycle: maturing)
    # Handle default column additions (eg, soma_joinid) and error checking on
    # required columns
    # @return An [`arrow::Schema`], which may be modified by the addition of
    # required columns.
    validate_schema = function(schema, index_column_names) {
      stopifnot(
        "'schema' must be a valid Arrow schema" = is_arrow_schema(schema),
        "'index_column_names' must be a non-empty character vector" = is.character(
          index_column_names
        ) &&
          length(index_column_names) > 0,
        "All 'index_column_names' must be defined in the 'schema'" = assert_subset(
          index_column_names,
          schema$names,
          "indexed field"
        ),
        "Column names must not start with reserved prefix 'soma_'" = all(
          !startsWith(setdiff(schema$names, "soma_joinid"), "soma_")
        ) ||
          isTRUE(getOption("tiledbsoma.write_soma.internal", default = FALSE))
      )

      # Add soma_joinid column if not present
      if ("soma_joinid" %in% schema$names) {
        stopifnot(
          "soma_joinid field must be of type Arrow int64" = schema$GetFieldByName(
            "soma_joinid"
          )$type ==
            arrow::int64()
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
      if (!.is_domain(new_domain, dimnames)) {
        stop(
          "new_domain must be a named list, with values being 2-element",
          " vectors or NULL, with names the same as the dataframe's",
          " index-column names"
        )
      }

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
        if (
          is.null(new_domain[[dimname]]) &&
            full_schema[[dimname]]$type$ToString() %in%
              c("string", "large_string", "utf8", "large_utf8")
        ) {
          ordered_new_domain[[dimname]] <- c("", "")
        } else {
          ordered_new_domain[[dimname]] <- new_domain[[dimname]]
        }
      }

      pyarrow_table <- arrow::arrow_table(
        as.data.frame(ordered_new_domain),
        schema = dim_schema
      )

      # We transfer to the arrow table via a pair of array and schema pointers
      dnaap <- nanoarrow::nanoarrow_allocate_array()
      dnasp <- nanoarrow::nanoarrow_allocate_schema()
      arrow::as_record_batch(pyarrow_table)$export_to_c(dnaap, dnasp)

      return(list(array = dnaap, schema = dnasp))
    }
  )
)
