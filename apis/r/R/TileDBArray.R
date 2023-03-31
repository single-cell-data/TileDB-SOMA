#' TileDB Array Base Class
#'
#' @description Base class for representing an individual TileDB array.
#' (lifecycle: experimental)
#'
#' @export
TileDBArray <- R6::R6Class(
  classname = "TileDBArray",
  inherit = TileDBObject,
  public = list(

    #' @description Print summary of the array. (lifecycle: experimental)
    print = function() {
      super$print()
      if (self$exists()) {
        cat("  dimensions:", string_collapse(self$dimnames()), "\n")
        cat("  attributes:", string_collapse(self$attrnames()), "\n")
      }
    },

    #' @description Return a [`TileDBArray`] object (lifecycle: experimental)
    #' @param ... Optional arguments to pass to `tiledb::tiledb_array()`
    #' @return A [`tiledb::tiledb_array`] object.
    tiledb_array = function(...) {
      args <- list(...)
      args$uri <- self$uri
      args$query_type <- "READ"
      args$query_layout <- "UNORDERED"
      args$ctx <- self$tiledbsoma_ctx$context()
      do.call(tiledb::tiledb_array, args)
    },

    #' @description Retrieve metadata from the TileDB array. (lifecycle: experimental)
    #' @param key The name of the metadata attribute to retrieve.
    #' @param prefix Filter metadata using an optional prefix. Ignored if `key`
    #'   is not NULL.
    #' @return A list of metadata values.
    get_metadata = function(key = NULL, prefix = NULL) {
      private$open("READ")
      on.exit(self$close())

      if (!is.null(key)) {
        metadata <- tiledb::tiledb_get_metadata(self$object, key)
      } else {
        # coerce tiledb_metadata to list
        metadata <- unclass(tiledb::tiledb_get_all_metadata(self$object))
        if (!is.null(prefix)) {
          metadata <- metadata[string_starts_with(names(metadata), prefix)]
        }
      }
      return(metadata)
    },

    #' @description Add list of metadata to the specified TileDB array. (lifecycle: experimental)
    #' @param metadata Named list of metadata to add.
    #' @return NULL
    set_metadata = function(metadata) {
      stopifnot(
        "Metadata must be a named list" = is_named_list(metadata)
      )
      private$open("WRITE")
      on.exit(self$close())

      dev_null <- mapply(
        FUN = tiledb::tiledb_put_metadata,
        key = names(metadata),
        val = metadata,
        MoreArgs = list(arr = self$object),
        SIMPLIFY = FALSE
      )
    },

    #' @description Retrieve the array schema as an Arrow schema (lifecycle: experimental)
    #' @return A [`arrow::schema`] object
    schema = function() {
      arrow_schema_from_tiledb_schema(tiledb::schema(self$object))
    },

    #' @description Retrieve the array schema as TileDB schema (lifecycle: experimental)
    #' @return A [`tiledb::tiledb_array_schema`] object
    tiledb_schema = function() {
      tiledb::schema(self$object)
    },

    #' @description Retrieve the array dimensions (lifecycle: experimental)
    #' @return A named list of [`tiledb::tiledb_dim`] objects
    dimensions = function() {
      dims <- tiledb::dimensions(self$tiledb_schema())
      setNames(dims, nm = vapply_char(dims, tiledb::name))
    },

    #' @description Retrieve the shape, i.e. the length of each dimension (lifecycle: experimental)
    #' @return A named vector of dimension length (and the same type as the dimension)
    shape = function() {
      as.integer64(shape(
        self$uri,
        config=as.character(tiledb::config(self$tiledbsoma_ctx$context()))
      ))
    },

    #' @description Retrieve number of dimensions (lifecycle: experimental)
    #' @return A scalar with the number of dimensions
    ndim = function() {
      dims <- tiledb::dimensions(self$tiledb_schema())
      length(dims)
    },

    #' @description Retrieve the array attributes (lifecycle: experimental)
    #' @return A list of [`tiledb::tiledb_attr`] objects
    attributes = function() {
      tiledb::attrs(self$tiledb_schema())
    },

    #' @description Retrieve dimension names (lifecycle: experimental)
    #' @return A character vector with the array's dimension names
    dimnames = function() {
      vapply(
        self$dimensions(),
        FUN = tiledb::name,
        FUN.VALUE = vector("character", 1L),
        USE.NAMES = FALSE
      )
    },

    #' @description Retrieve attribute names (lifecycle: experimental)
    #' @return A character vector with the array's attribute names
    attrnames = function() {
      vapply(
        self$attributes(),
        FUN = tiledb::name,
        FUN.VALUE = vector("character", 1L),
        USE.NAMES = FALSE
      )
    },

    #' @description Retrieve the names of all columns, including dimensions and
    #' attributes (lifecycle: experimental)
    #' @return A character vector with the array's column names
    colnames = function() {
      c(self$dimnames(), self$attrnames())
    },

    #' @description Retrieve names of index (dimension) columns (lifecycle: experimental)
    #' @return A character vector with the array index (dimension) names
    index_column_names = function() {
      self$dimnames()
    },

    #' @description Get number of fragments in the array (lifecycle: experimental)
    fragment_count = function() {
      tiledb::tiledb_fragment_info_get_num(
        tiledb::tiledb_fragment_info(self$uri)
      )
    },

    #' @description Set dimension values to slice from the array. (lifecycle: experimental)
    #' @param dims a named list of character vectors. Each name must correspond
    #' to an array dimension. The character vectors within each element are used
    #' to set the arrays selected ranges for each corresponding dimension.
    #' @param attr_filter a TileDB query condition for attribute filtering.
    #' pushdown.
    set_query = function(dims = NULL, attr_filter = NULL) {

      if (!is.null(dims)) {
        stopifnot(
          "'dims' must be a named list of character vectors" =
            is_named_list(dims) && all(vapply_lgl(dims, is_character_or_null)),
          assert_subset(names(dims), self$dimnames(), type = "dimension")
        )

        # list of dimensions to slice discarding NULL elements
        dims <- modifyList(list(), dims)

        if (!is_empty(dims)) {
          # Convert each dim vector to a two-column matrix where each row
          # describes one pair of minimum and maximum values.
          tiledb::selected_ranges(private$tiledb_object) <- lapply(
            X = dims,
            FUN = function(x) unname(cbind(x, x))
          )
        }
      }

      # We have to add special handling for attr_filter to cover the case where
      # 1) TileDBArray$set_query() was called directly and attr_filter is an
      # unevaluated expression, and 2) when TileDBArray$set_query() was called
      # indirectly (via AnnotationGroup$set_query()) and attr_filter has been
      # captured and converted to a a character vector.

      # capture error thrown if attr_filter is an unevaluated expression
      # suppress warning: restarting interrupted promise evaluation
      is_character_expression <- suppressWarnings(
          try(is.character(attr_filter), silent = TRUE)
      )

      if (inherits(is_character_expression, "try-error")) {
          # attr_filter is an unevaluated expression
          captured_filter <- substitute(attr_filter)
      } else if (is_character_expression) {
          # attr_filter has already been captured and converted to a char vector
          if (attr_filter == "NULL") return(NULL)
          captured_filter <- str2lang(attr_filter)
      } else if (is.null(attr_filter)) {
          return(NULL)
      } else {
          stop("'attr_filter' is not a valid expression", call. = FALSE)
      }

      tiledb::query_condition(private$tiledb_object) <- do.call(
        what = tiledb::parse_query_condition,
        args = list(expr = captured_filter, ta = self$object)
      )
    },

    #' @description Reset the query. By default both dimension ranges and
    #' attribute filters are cleared. (lifecycle: experimental)
    #' @param dims Clear the defined dimension ranges?
    #' @param attr_filter Clear the defined attribute filters?
    #' @return NULL
    reset_query = function(dims = TRUE, attr_filter = TRUE) {
      stopifnot(
        "'dims' must be a logical" = is.logical(dims),
        "'attr_filter' must be a logical" = is.logical(attr_filter),
        "Nothing will be reset if 'dims' and 'attr_filter' are both FALSE"
          = isTRUE(dims) || isTRUE(attr_filter)
      )
      if (isTRUE(dims)) {
        tiledb::selected_ranges(private$tiledb_object) <- list()
      }
      if (isTRUE(attr_filter)) {
        tiledb::query_condition(private$tiledb_object) <- new(
          "tiledb_query_condition"
        )
      }
    },

    #' @description Close the SOMA object.
    #' @return The object, invisibly
    close = function() {
      spdl::debug("Closing {} '{}'", self$class(), self$uri)
      invisible(tiledb::tiledb_array_close(self$object))
    }

  ),

  private = list(

    # Once the array has been created this initializes the TileDB array object
    # and stores the reference in private$tiledb_object.
    initialize_object = function() {
      private$tiledb_object <- tiledb::tiledb_array(
        uri = self$uri,
        ctx = self$tiledbsoma_ctx$context(),
        query_layout = "UNORDERED"
      )
      self$close()
    },

    open = function(mode) {
      mode <- match.arg(mode, c("READ", "WRITE"))
      spdl::debug(
        "Opening {} '{}' in {} mode", self$class(), self$uri, mode
      )
      invisible(tiledb::tiledb_array_open(self$object, type = mode))
    }
  )
)
