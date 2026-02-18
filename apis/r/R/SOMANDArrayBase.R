#' SOMA ND-Array Base Class
#'
#' @description Virtual base class to add ND-array-specific functionality to the
#' \code{\link{SOMAArrayBase}} class (lifecycle: maturing).
#'
#' @param check_only If true, does not apply the operation, but only reports
#' whether it would have succeeded.
#'
#' @keywords internal
#'
#' @export
#'
#' @seealso Derived classes: \code{\link{SOMADenseNDArray}},
#' \code{\link{SOMASparseNDArray}}
#'
SOMANDArrayBase <- R6::R6Class(
  classname = "SOMANDArrayBase",
  inherit = SOMAArrayBase,
  public = list(
    #' @description Create a SOMA NDArray named with the URI
    #' (lifecycle: maturing).\cr
    #' \cr
    #' \strong{Note}: \code{$create()} is considered internal and should not be
    #' called directly; use factory functions
    #' (eg. \code{\link{SOMASparseNDArrayCreate}()}) instead.
    #'
    #' @param type An \link[arrow:data-type]{Arrow type} defining the type
    #' of each element in the array.
    #' @param shape a vector of integers defining the shape of the array.
    #' @template param-platform-config
    #'
    #' @return Returns \code{self}.
    #'
    create = function(type, shape, platform_config = NULL) {
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

      ## .is_sparse field is being set by dense and sparse private initialisers, respectively
      private$.type <- type # Arrow schema type of data

      dom_ext_tbl <- get_domain_and_extent_array(shape, self$is_sparse())

      # Parse the tiledb/create subkeys of the platform_config into a handy,
      # typed, queryable data structure.
      tiledb_create_options <- TileDBCreateOptions$new(platform_config)

      ## create array
      # ctxptr <- self$tiledbsoma_ctx$context()
      sparse <- if (inherits(self, "SOMASparseNDArray")) {
        TRUE
      } else if (inherits(self, "SOMADenseNDArray")) {
        FALSE
      } else {
        stop("Unknown SOMA array type: ", self$class(), call. = FALSE)
      }
      createSchemaForNDArray(
        uri = self$uri,
        format = as_nanoarrow_schema(arrow::schema(arrow::field("soma_data", type)))$children$soma_data$format,
        shape = as.integer64(shape),
        soma_type = if (sparse) "SOMASparseNDArray" else "SOMADenseNDArray",
        pclst = tiledb_create_options$to_list(FALSE),
        ctxxp = private$.context$handle,
        tsvec = self$.tiledb_timestamp_range
      )
      # private$write_object_type_metadata(timestamps)  ## FIXME: temp. commented out -- can this be removed overall?

      self$reopen("WRITE", tiledb_timestamp = self$tiledb_timestamp)
      return(self)
    },

    ## needed eg after open() to set (Arrow) type
    #' @description Sets a cache value for the datatype.
    #'
    #' @section Lifecycle:
    #' As of \pkg{tiledbsoma} 2.1.0, \code{$set_data_type()} is deprecated; this
    #' functionality is no longer needed as \code{libtiledbsoma} now accurately
    #' sets the TileDB data type upon array opening
    #'
    #' @param type A character value describing the TileDB data type.
    #'
    set_data_type = function(type) {
      .deprecate(
        when = "2.1.0",
        what = sprintf("%s$set_data_type()", class(self)[1L])
      )
      soma_debug(sprintf(
        "[SOMANDArrayBase::set_data_type] caching type %s",
        type$ToString()
      ))
      private$.type <- type
    },

    #' @description Test if the array has the upgraded resizeable shape feature
    #' from TileDB-SOMA 1.15, the array was created with this support, or it has
    #' had \code{$upgrade_domain()} applied to it (lifecycle: maturing).
    #'
    #' @return Returns \code{TRUE} if the array has the upgraded resizable
    #' shape feature; otherwise, returns \code{FALSE}.
    #'
    tiledbsoma_has_upgraded_shape = function() {
      has_current_domain(self$uri, private$.context$handle)
    },

    #' @description Increases the shape of the array as specified, up to the hard
    #' limit which is \code{maxshape}. Raises an error if the new shape is less
    #' than the current shape or exceeds \code{maxshape} in any dimension. Also
    #' raises an error if the array doesn't already have a shape; in that case
    #' please call \code{$tiledbsoma_upgrade_shape()} (lifecycle: maturing).
    #' @param new_shape An integerish vector of the same length as the array's
    #' \code{$ndim()}.
    #'
    #' @return If \code{check_only}, returns the empty string if no error is
    #' detected, else a description of the error. Otherwise, invisibly returns
    #' \code{NULL}.
    #'
    resize = function(new_shape, check_only = FALSE) {
      stopifnot(
        "'new_shape' must be a vector of integerish values, of the same length as maxshape" = rlang::is_integerish(
          new_shape,
          n = self$ndim()
        ) ||
          (bit64::is.integer64(new_shape) && length(new_shape) == self$ndim())
      )
      # Checking slotwise new shape >= old shape, and <= max_shape, is already done in libtiledbsoma

      reason_string <- resize(
        uri = self$uri,
        new_shape = new_shape,
        function_name_for_messages = .name_of_function(),
        check_only = check_only,
        ctxxp = private$.context$handle
      )

      if (isTRUE(check_only)) {
        return(reason_string)
      }

      # The return value from resize without check_only is always "", or it
      # raises an error trying.
      return(invisible(NULL))
    },

    #' @description Allows the array to have a resizeable shape as described in
    #' the TileDB-SOMA 1.15 release notes. Raises an error if the shape exceeds
    #' \code{maxshape} in any dimension, or if the array already has a shape.
    #' The methods \code{$tiledbsoma_upgrade_shape()} and \code{$resize()} are
    #' very similar: the former must be called on a pre-1.15 array the first
    #' time a shape is set on it; the latter must be used for subsequent resizes
    #' on any array which already has upgraded shape (lifecycle: maturing).
    #'
    #' @param shape An integerish vector of the same length as the array's
    #' \code{$ndim()}.
    #'
    #' @return If \code{check_only}, returns the empty string if no error is
    #' detected, else a description of the error. Otherwise, invisibly returns
    #' \code{NULL}.
    #'
    tiledbsoma_upgrade_shape = function(shape, check_only = FALSE) {
      stopifnot(
        "'shape' must be a vector of integerish values, of the same length as maxshape" = rlang::is_integerish(
          shape,
          n = self$ndim()
        ) ||
          (bit64::is.integer64(shape) && length(shape) == self$ndim())
      )
      # Checking slotwise new shape >= old shape, and <= max_shape, is already done in libtiledbsoma

      reason_string <- tiledbsoma_upgrade_shape(
        uri = self$uri,
        new_shape = shape,
        function_name_for_messages = .name_of_function(),
        check_only = check_only,
        ctxxp = private$.context$handle
      )
      if (isTRUE(check_only)) {
        return(reason_string)
      }

      # The return value from tiledbsoma_upgrade_shape without check_only is
      # always "", or it raises an error trying.
      return(invisible(NULL))
    }
  ),
  private = list(
    .type = NULL,

    # Given a user-specified shape along a particular dimension, returns a named
    # list containing name, capacity, and extent elements.
    .dim_capacity_and_extent = function(name, shape = NULL, create_options) {
      stop(
        ".dim_capacity_and_extent() must be implemented by child class",
        call. = FALSE
      )
    },

    #  @description Converts a list of vectors corresponding to coords to a
    #  format acceptable for mq_setup and soma_array_reader
    .convert_coords = function(coords) {
      # Ensure coords is a named list, use to select dim points
      stopifnot(
        "'coords' must be a list" = is.list(coords) && length(coords),
        "'coords' must be a list integerish vectors" = all(vapply(
          X = coords,
          FUN = function(x) {
            if (is.null(x)) {
              return(TRUE)
            }
            return(
              (is.null(dim(x)) && !is.factor(x)) &&
                (rlang::is_integerish(x, finite = TRUE) ||
                  (bit64::is.integer64(x) && all(is.finite(x)))) &&
                length(x) &&
                all(x >= 0L)
            )
          },
          FUN.VALUE = logical(length = 1L),
          USE.NAMES = FALSE
        )),
        "'coords' if unnamed must have length of dim names, else if named names must match dim names" = ifelse(
          test = is.null(names(coords)),
          yes = length(coords) == length(self$dimnames()),
          no = all(names(coords) %in% self$dimnames())
        )
      )

      # Remove NULL-entries from coords
      coords <- Filter(Negate(is.null), coords)
      if (!length(coords)) {
        return(NULL)
      }

      # If unnamed, set names
      if (is.null(names(coords))) {
        names(coords) <- self$dimnames()
      }

      # Convert to integer64 to match dimension type
      return(sapply(
        coords,
        FUN = bit64::as.integer64,
        simplify = FALSE,
        USE.NAMES = TRUE
      ))
    },

    #  @description Converts a vector of ints into a vector of int64 in a format
    #  acceptable for libtiledbsoma

    .convert_shape_argument = function(new_shape) {
      # ensure new_shape is an integerish vector
      stopifnot(
        "'new_shape' must be an integerish vector with the same length as the array's maxshape" = rlang::is_integerish(
          new_shape,
          n = self$ndim(),
          finite = TRUE
        ) ||
          (bit64::is.integer64(new_shape) &&
            length(new_shape) == self$ndim() &&
            all(is.finite(new_shape)))
      )

      # convert integer to integer64 to match dimension type
      return(bit64::as.integer64(new_shape))
    }
  )
)
