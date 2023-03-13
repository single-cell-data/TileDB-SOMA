#' TileDB Array Base Class
#'
#' @description
#' Base class for representing an individual TileDB array.
#'
#' @details
#' ## Initialization
#' Initializing a `TileDBArray` object does not automatically create a new array
#' at the specified `uri` if one does not already exist because we don't know
#' what the schema will be. Arrays are only created by child classes, which
#' populate the private `create_empty_array()` and `ingest_data()` methods.
#' @export
TileDBArray <- R6::R6Class(
  classname = "TileDBArray",
  inherit = TileDBObject,
  public = list(

    #' @description Create a new TileDBArray object.
    #' @param uri URI for the TileDB array
    #' @param verbose Print status messages
    #' @param config optional configuration
    #' @param ctx optional tiledb context
    initialize = function(uri, verbose = TRUE, config = NULL, ctx = NULL) {
      super$initialize(uri, verbose, config, ctx)

      if (self$exists()) {
        msg <- sprintf("Found existing %s at '%s'", self$class(), self$uri)
        private$initialize_object()
        # Check for legacy validity mode metadata tag
        toggle_tiledb_legacy_mode_if_needed(self$object, self$verbose)
      } else {
        msg <- sprintf("No %s found at '%s'", self$class(), self$uri)
      }
      if (self$verbose) message(msg)
      return(self)
    },

    #' @description Print summary of the array.
    print = function() {
      super$print()
      if (self$exists()) {
        cat("  dimensions:", string_collapse(self$dimnames()), "\n")
        cat("  attributes:", string_collapse(self$attrnames()), "\n")
      }
    },

    #' @description Check if the array exists.
    #' @return TRUE if the array exists, FALSE otherwise.
    array_exists = function() {
      .Deprecated(
        new = "exists()",
        old = "array_exists()"
      )
      self$exists()
    },

    #' @description Return a [`TileDBArray`] object
    #' @param ... Optional arguments to pass to `tiledb::tiledb_array()`
    #' @return A [`tiledb::tiledb_array`] object.
    tiledb_array = function(...) {
      args <- list(...)
      args$uri <- self$uri
      args$query_type <- "READ"
      args$query_layout <- "UNORDERED"
      args$ctx <- self$ctx
      do.call(tiledb::tiledb_array, args)
    },

    #' @description Retrieve metadata from the TileDB array.
    #' @param key The name of the metadata attribute to retrieve.
    #' @param prefix Filter metadata using an optional prefix. Ignored if `key`
    #'   is not NULL.
    #' @return A list of metadata values.
    get_metadata = function(key = NULL, prefix = NULL) {
      on.exit(private$close())
      private$open("READ")
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

    #' @description Add list of metadata to the specified TileDB array.
    #' @param metadata Named list of metadata to add.
    #' @param prefix Optional prefix to add to the metadata attribute names.
    #' @return NULL
    add_metadata = function(metadata, prefix = "") {
      stopifnot(
        "Metadata must be a named list" = is_named_list(metadata)
      )
      on.exit(private$close())
      if (self$verbose) {
        message(sprintf("Adding %i metadata keys to array", length(metadata)))
      }
      private$open("WRITE")
      mapply(
        FUN = tiledb::tiledb_put_metadata,
        key = paste0(prefix, names(metadata)),
        val = metadata,
        MoreArgs = list(arr = self$object),
        SIMPLIFY = FALSE
      )
    },

    #' @description Retrieve the array schema
    #' @return A [`tiledb::tiledb_array_schema`] object
    schema = function() {
      tiledb::schema(self$object)
    },

    #' @description Retrieve the array dimensions
    #' @return A list of [`tiledb::tiledb_dim`] objects
    dimensions = function() {
      tiledb::dimensions(self$schema())
    },

    #' @description Retrieve the array attributes
    #' @return A list of [`tiledb::tiledb_attr`] objects
    attributes = function() {
      tiledb::attrs(self$schema())
    },

    #' @description Retrieve dimension names
    #' @return A character vector with the array's dimension names
    dimnames = function() {
      vapply(
        self$dimensions(),
        FUN = tiledb::name,
        FUN.VALUE = vector("character", 1L)
      )
    },

    #' @description Get number of fragments in the array
    fragment_count = function() {
      tiledb::tiledb_fragment_info_get_num(
        tiledb::tiledb_fragment_info(self$uri)
      )
    },

    #' @description Retrieve attribute names
    #' @return A character vector with the array's attribute names
    attrnames = function() {
      vapply(
        self$attributes(),
        FUN = tiledb::name,
        FUN.VALUE = vector("character", 1L),
        USE.NAMES = FALSE
      )
    },

    #' @description Set dimension values to slice from the array.
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
          stop("'attr_filter' is not a valid expression")
      }

      tiledb::query_condition(private$tiledb_object) <- do.call(
        what = tiledb::parse_query_condition,
        args = list(expr = captured_filter, ta = self$object)
      )
    },

    #' @description Reset the query. By default both dimension ranges and
    #' attribute filters are cleared.
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
    }
  ),

  private = list(

    # Once the array has been created this initializes the TileDB array object
    # and stores the reference in private$tiledb_object.
    initialize_object = function() {
      private$tiledb_object <- tiledb::tiledb_array(
        uri = self$uri,
        ctx = self$ctx,
        query_layout = "UNORDERED"
      )

      private$close()
    },

    write_object_type_metadata = function() {
      meta <- list()
      meta[[SOMA_OBJECT_TYPE_METADATA_KEY]] <- class(self)[1]
      meta[[SOMA_ENCODING_VERSION_METADATA_KEY]] <- SOMA_ENCODING_VERSION
      meta[[SOMA_LEGACY_VALIDITY_KEY]] <- SOMA_LEGACY_VALIDITY
      self$add_metadata(meta) # TileDBArray or TileDBGroup
    },

    # @description Create empty TileDB array.
    create_empty_array = function() return(NULL),

    open = function(mode) {
      mode <- match.arg(mode, c("READ", "WRITE"))
      invisible(tiledb::tiledb_array_open(self$object, type = mode))
    },

    close = function() {
      invisible(tiledb::tiledb_array_close(self$object))
    },

    # @description Ingest data into the TileDB array.
    ingest_data = function() return(NULL),

    # @description Retrieve data from the TileDB array
    # @param batch_mode logical, if `TRUE`, batch query mode is enabled, which
    # provides the ability to detect partial query results and resubmit until
    # all results are retrieved.
    # @param return_as Data can be read in as a `list` (default), `array`,
    # `matrix`, `data.frame`, `data.table` or `tibble`.
    read_data = function(attrs = NULL, batch_mode = FALSE, return_as = NULL) {
      if (self$verbose) {
        message(
          sprintf("Reading %s into memory from '%s'", self$class(), self$uri)
        )
      }
      arr <- self$object
      tiledb::attrs(arr) <- attrs %||% character()
      tiledb::return_as(arr) <- return_as %||% "asis"

      if (batch_mode) {
        if (self$verbose) message("...reading in batches")
        batcher <- tiledb::createBatched(arr)
        results <- list()
        i <- 1
        while (isFALSE(tiledb::completedBatched(batcher))) {
          if (self$verbose) message(sprintf("...retrieving batch %d", i))
          results[[i]] <- tiledb::fetchBatched(arr, batcher)
          i <- i + 1
        }

        # TODO: currently tiledb-r's batched reader ignores return_as and a
        # data.frame is always returned. When this is addressed we'll need to
        # add class-specific concatenation logic here.
        results <- vctrs::vec_rbind(!!!results)
      } else {
        results <- arr[]
      }
      results
    }
  )
)
