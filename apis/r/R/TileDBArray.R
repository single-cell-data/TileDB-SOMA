#' TileDB Array Base Class
#'
#' @description Virtual base class for representing an individual TileDB array
#' (lifecycle: maturing).
#'
#' @keywords internal
#'
#' @export
#'
#' @seealso Derived classes: \code{\link{SOMAArray}}
#'
TileDBArray <- R6::R6Class(
  classname = "TileDBArray",
  inherit = TileDBObject,
  public = list(

    #' @description Open the SOMA object for read or write.
    #'
    #' @param mode Mode to open in; defaults to \code{READ}.
    #' @param internal_use_only Character value to signal this is a 'permitted'
    #' call, as \code{open()} is considered internal and should not be called
    #' directly.
    #'
    #' @return Returns \code{self}
    #'
    open = function(mode = c("READ", "WRITE"), internal_use_only = NULL) {
      mode <- match.arg(mode)
      if (is.null(internal_use_only) || internal_use_only != "allowed_use") {
        stop(
          "Use of the open() method is for internal use only.\n",
          "  - Consider using a factory method as e.g. 'SOMADataFrameOpen()'.",
          call. = FALSE
        )
      }

      private$.mode <- mode
      if (is.null(self$tiledb_timestamp)) {
        spdl::debug("[TileDBArray$open] Opening {} '{}' in {} mode", self$class(), self$uri, mode)
        private$.tiledb_array <- tiledb::tiledb_array_open(self$object, type = mode)
      } else {
        if (is.null(internal_use_only)) stopifnot("tiledb_timestamp not yet supported for WRITE mode" = mode == "READ")
        spdl::debug(
          "[TileDBArray$open] Opening {} '{}' in {} mode at ({})",
          self$class(),
          self$uri,
          mode,
          self$tiledb_timestamp %||% "now"
        )
        # private$.tiledb_array <- tiledb::tiledb_array_open_at(self$object, type = mode,
        #                                                      timestamp = self$tiledb_timestamp)
      }

      ## TODO -- cannot do here while needed for array case does not work for data frame case
      # private$.type <- arrow_type_from_tiledb_type(tdbtype)

      private$update_metadata_cache()
      self
    },

    #' @description Close the SOMA object.
    #'
    #' @return Invisibly returns \code{self}.
    #'
    close = function() {
      spdl::debug("[TileDBArray$close] Closing {} '{}'", self$class(), self$uri)
      private$.mode <- "CLOSED"
      tiledb::tiledb_array_close(self$object)
      invisible(self)
    },

    #' @description Print summary of the array (lifecycle: maturing).
    #'
    print = function() {
      super$print()
      if (self$exists()) {
        cat("  dimensions:", string_collapse(self$dimnames()), "\n")
        cat("  attributes:", string_collapse(self$attrnames()), "\n")
      }
    },

    #' @description Return a \code{\link[tiledb]{tiledb_array}} object
    #' (lifecycle: maturing).
    #'
    #' @param ... Optional arguments to pass to
    #' \code{\link[tiledb:tiledb_array]{tiledb::tiledb_array}()}.
    #'
    #' @return A \code{\link[tiledb]{tiledb_array}} object.
    #'
    tiledb_array = function(...) {
      args <- list(...)
      args$uri <- self$uri
      args$query_type <- self$.mode
      args$query_layout <- "UNORDERED"
      args$ctx <- self$tiledbsoma_ctx$context()
      spdl::debug("[TileDBArray$tiledb_array] ctor uri {} mode {} layout {}", args$uri, args$query_type, args$query_layout)
      do.call(tiledb::tiledb_array, args)
    },

    #' @description Retrieve metadata from the TileDB array
    #' (lifecycle: maturing).
    #'
    #' @param key The name of the metadata attribute to retrieve.
    #'
    #' @return A list of metadata values.
    #'
    get_metadata = function(key = NULL) {
      private$check_open_for_read_or_write()

      spdl::debug("[TileDBArray$get_metadata] Retrieving metadata for {} '{}'", self$class(), self$uri)
      private$fill_metadata_cache_if_null()
      if (!is.null(key)) {
        val <- private$.metadata_cache[[key]]
        if (is.list(val)) val <- unlist(val)
        val
      } else {
        private$.metadata_cache
      }
    },

    #' @description Add list of metadata to the specified TileDB array
    #' (lifecycle: maturing).
    #'
    #' @param metadata Named list of metadata to add.
    #'
    #' @return Invisibly returns \code{NULL}.
    #'
    set_metadata = function(metadata) {
      stopifnot(
        "Metadata must be a named list" = is_named_list(metadata)
      )

      # private$check_open_for_write()

      for (nm in names(metadata)) {
        val <- metadata[[nm]]
        spdl::debug("[TileDBArray$set_metadata] setting key {} to {} ({})", nm, val, class(val))
        set_metadata(
          uri = self$uri,
          key = nm,
          valuesxp = val,
          type = class(val),
          is_array = TRUE,
          ctxxp = soma_context(),
          tsvec = self$.tiledb_timestamp_range
        )
      }

      dev_null <- mapply(
        FUN = private$add_cached_metadata,
        key = names(metadata),
        val = metadata,
        SIMPLIFY = FALSE
      )
    },

    #' @description Retrieve the array schema as an Arrow schema
    #' (lifecycle: maturing).
    #'
    #' @return An \link[arrow:Schema]{Arrow schema} object.
    #'
    schema = function() {
      return(arrow::as_schema(c_schema(self$uri, private$.soma_context)))
    },

    #' @description Retrieve the array schema as TileDB schema
    #' (lifecycle: maturing).
    #'
    #' @return A \code{\link[tiledb]{tiledb_array_schema}} object.
    #'
    tiledb_schema = function() {
      tiledb::schema(self$object)
    },

    #' @description Retrieve the array dimensions (lifecycle: maturing).
    #'
    #' @return A named list of \code{\link[tiledb]{tiledb_dim}} objects.
    #'
    dimensions = function() {
      dims <- tiledb::dimensions(self$tiledb_schema())
      return(stats::setNames(dims, nm = vapply_char(dims, tiledb::name)))
    },

    #' @description Retrieve the shape, i.e. the capacity of each dimension.
    #' Attempted reads and writes outside the \code{shape} will result in a
    #' runtime error: this is the purpose of \code{shape}. This will not
    #' necessarily match the bounds of occupied cells within the array.
    #' Using \code{$resize()}, this may be increased up to the hard limit which
    #' \code{$maxshape()} reports (lifecycle: maturing).
    #'
    #' @return A named vector of dimension length
    #' (and the same type as the dimension).
    #'
    shape = function() {
      return(bit64::as.integer64(shape(self$uri, private$.soma_context)))
    },

    #' @description Retrieve the hard limit up to which the array may be resized
    #' using the \code{$resize()} method (lifecycle: maturing).
    #'
    #' @return A named vector of dimension length
    #' (and the same type as the dimension).
    #'
    maxshape = function() {
      return(bit64::as.integer64(maxshape(self$uri, private$.soma_context)))
    },

    #' @description Returns a named list of minimum/maximum pairs, one per
    #' index column, which are the smallest and largest values written on that
    #' index column.
    #'
    #' @param index1 Return the non-empty domain with 1-based indices.
    #' @param max_only Return only the max value per dimension, and return
    #' this as a vector. Names are dropped (lifecycle: maturing).
    #'
    #' @return Named list of minimum/maximum values, or integer vector
    #' of maximum values.
    #'
    non_empty_domain = function(index1 = FALSE, max_only = FALSE) {
      retval <- as.list(
        arrow::as_record_batch(
          arrow::as_arrow_table(
            non_empty_domain(self$uri, private$.soma_context)
          )
        )
      )
      if (index1) {
        retval <- lapply(retval, function(c) {
          c + 1
        })
      }
      if (max_only) {
        # No vapply options since SOMADataFrame can have varying types.
        retval <- unname(unlist(lapply(retval, function(e) {
          e[[2]]
        })))
      }
      return(retval)
    },

    #' @description Retrieve number of dimensions (lifecycle: maturing).
    #'
    #' @return A scalar with the number of dimensions.
    #'
    ndim = function() {
      ndim(self$uri, private$.soma_context)
    },

    #' @description Retrieve the array attributes (lifecycle: maturing).
    #'
    #' @return A list of \code{\link[tiledb]{tiledb_attr}} objects.
    #'
    attributes = function() {
      tiledb::attrs(self$tiledb_schema())
    },

    #' @description Retrieve dimension names (lifecycle: maturing).
    #'
    #' @return A character vector with the array's dimension names.
    #'
    dimnames = function() {
      c_dimnames(self$uri, private$.soma_context)
    },

    #' @description Retrieve attribute names (lifecycle: maturing).
    #'
    #' @return A character vector with the array's attribute names.
    #'
    attrnames = function() {
      c_attrnames(self$uri, private$.soma_context)
    },

    #' @description Retrieve the names of all columns, including dimensions and
    #' attributes (lifecycle: maturing).
    #'
    #' @return A character vector with the array's column names.
    #'
    colnames = function() {
      c(self$dimnames(), self$attrnames())
    },

    #' @description Retrieve names of index (dimension) columns
    #' (lifecycle: maturing).
    #'
    #' @return A character vector with the array index (dimension) names.
    #'
    index_column_names = function() {
      self$dimnames()
    }
  ),
  active = list(
    #' @field object Access the underlying TileB object directly (either a
    #' \code{\link[tiledb]{tiledb_array}} or \code{\link[tiledb]{tiledb_group}}).
    #'
    object = function(value) {
      if (!missing(value)) {
        stop(sprintf("'%s' is a read-only field.", "object"), call. = FALSE)
      }
      # If the array was created after the object was instantiated, we need to
      # initialize private$.tiledb_array
      if (is.null(private$.tiledb_array)) {
        private$initialize_object()
      }
      private$.tiledb_array
    }
  ),
  private = list(

    # Internal pointer to the TileDB array.
    #
    # Important implementation note:
    # * In TileDB-R there is an unopened handle obtained by tiledb::tiledb_array, which takes
    #   a URI as its argument.
    # * One may then open and close this using tiledb::tiledb_array_open (for read or write)
    #   and tiledb::tiledb_array_close, which take a tiledb_array handle as their first argument.
    #
    # However, for groups:
    # * tiledb::tiledb_group and tiledb::group_open both return an object opened for read or write.
    # * Therefore for groups we cannot imitate the behavior for arrays.
    #
    # For this reason there is a limit to how much handle-abstraction we can do in the TileDBObject
    # parent class. In particular, we cannot have a single .tiledb_object shared by both TileDBArray
    # and TileDBGroup.
    .tiledb_array = NULL,

    # Initially NULL, once the array is created or opened, this is populated
    # with a list that's empty or contains the array metadata. Since the SOMA
    # spec requires that we allow readback of array metadata even when the array
    # is open for write, but the TileDB layer underneath us does not, we must
    # have this cache.
    .metadata_cache = NULL,

    # Once the array has been created this initializes the TileDB array object
    # and stores the reference in private$.tiledb_array.
    initialize_object = function() {
      private$.tiledb_array <- tiledb::tiledb_array(
        uri = self$uri,
        ctx = self$tiledbsoma_ctx$context(),
        query_layout = "UNORDERED"
      )
    },

    # ----------------------------------------------------------------
    # Metadata-caching

    fill_metadata_cache_if_null = function() {
      if (is.null(private$.metadata_cache)) {
        private$update_metadata_cache()
      }
    },
    update_metadata_cache = function() {
      spdl::debug("[TileDBArray$update_metadata_cache] updating metadata cache for {} '{}' in {}", self$class(), self$uri, private$.mode)

      # See notes above -- at the TileDB implementation level, we cannot read array metadata
      # while the array is open for read, but at the SOMA application level we must support
      # this. Therefore if the array is opened for write and there is no cache populated then
      # we must open a temporary handle for read, to fill the cache.
      # array_handle <- private$.tiledb_array
      # if (private$.mode == "WRITE") {
      #  spdl::debug("[TileDBArray::update_metadata_cache] getting object")
      #  array_object <- tiledb::tiledb_array(self$uri, ctx = private$.tiledb_ctx)
      #  array_handle <- tiledb::tiledb_array_open(array_object, type = "READ")
      # }

      # if (isFALSE(tiledb::tiledb_array_is_open(array_handle))) {
      #  spdl::debug("[TileDBArray::update_metadata_cache] reopening object")
      #  array_handle <- tiledb::tiledb_array_open(array_handle, type = "READ")
      # }

      private$.metadata_cache <- get_all_metadata(self$uri, TRUE, soma_context())
      # print(str(private$.metadata_cache))
      # if (private$.mode == "WRITE") {
      #  tiledb::tiledb_array_close(array_handle)
      # }

      invisible(NULL)
    },
    add_cached_metadata = function(key, value) {
      if (is.null(private$.metadata_cache)) {
        private$.metadata_cache <- list()
      }
      private$.metadata_cache[[key]] <- value
    }
  )
)
