#' SOMA NDArray Base Class
#'
#' @description
#' Adds NDArray-specific functionality to the [`SOMAArrayBase`] class.
#' (lifecycle: maturing)
#'
#' @keywords internal
#' @export
#' @importFrom bit64 as.integer64

SOMANDArrayBase <- R6::R6Class(
  classname = "SOMANDArrayBase",
  inherit = SOMAArrayBase,
  public = list(

    #' @description Create a SOMA NDArray named with the URI. (lifecycle:
    #' experimental)
    #' @param type an [Arrow type][arrow::data-type] defining the type of each
    #' element in the array.
    #' @param shape a vector of integers defining the shape of the array.
    #' @template param-platform-config
    #' @param internal_use_only Character value to signal this is a 'permitted'
    #' call, as `create()` is considered internal and should not be called
    #' directly.
    create = function(
      type,
      shape,
      platform_config = NULL,
      internal_use_only = NULL
    ) {
      if (is.null(internal_use_only) || internal_use_only != "allowed_use") {
        stop(paste(
          "Use of the create() method is for internal use only. Consider using a",
          "factory method as e.g. 'SOMASparseNDArrayCreate()'."
        ), call. = FALSE)
      }

      ## .is_sparse field is being set by dense and sparse private initialisers, respectively
      private$.type <- type # Arrow schema type of data

      dom_ext_tbl <- get_domain_and_extent_array(shape, private$.is_sparse)

      # Parse the tiledb/create/ subkeys of the platform_config into a handy,
      # typed, queryable data structure.
      tiledb_create_options <- TileDBCreateOptions$new(platform_config)

      ## we transfer to the arrow table via a pair of array and schema pointers
      dnaap <- nanoarrow::nanoarrow_allocate_array()
      dnasp <- nanoarrow::nanoarrow_allocate_schema()
      arrow::as_record_batch(dom_ext_tbl)$export_to_c(dnaap, dnasp)

      ## we need a schema pointer to transfer the schema information
      ## so we first embed the (single column) 'type' into a schema and
      ## combine it with domain schema
      schema <- arrow::unify_schemas(
        arrow::schema(dom_ext_tbl),
        arrow::schema(arrow::field("soma_data", type))
      )
      nasp <- nanoarrow::nanoarrow_allocate_schema()
      schema$export_to_c(nasp)

      ## create array
      ctxptr <- super$tiledbsoma_ctx$context()
      createSchemaFromArrow(
        uri = self$uri,
        nasp = nasp,
        nadimap = dnaap,
        nadimsp = dnasp,
        sparse = private$.is_sparse,
        datatype = if (private$.is_sparse) "SOMASparseNDArray" else "SOMADenseNDArray",
        pclst = tiledb_create_options$to_list(FALSE),
        ctxxp = private$.soma_context,
        tsvec = self$.tiledb_timestamp_range
      )
      # private$write_object_type_metadata(timestamps)  ## FIXME: temp. commented out -- can this be removed overall?

      self$open("WRITE", internal_use_only = "allowed_use")
      self
    },

    ## needed eg after open() to set (Arrow) type
    #' @description Sets a cache value for the datatype (lifecycle: experimental)
    #' @param type A character value describing the TileDB data type
    set_data_type = function(type) {
      spdl::debug("[SOMANDArrayBase::set_data_type] caching type {}", type$ToString())
      private$.type <- type
    },

    #' @description Returns TRUE if the array has the upgraded resizeable shape
    #' feature from TileDB-SOMA 1.15: the array was created with this support,
    #' or it has had ``upgrade_domain`` applied to it.
    #' (lifecycle: maturing)
    #' @return Logical
    tiledbsoma_has_upgraded_shape = function() {
      has_current_domain(self$uri, private$.soma_context)
    },

    #' @description Increases the shape of the array as specified, up to the hard
    #' limit which is `maxshape`. Raises an error if the new shape is less than
    #' the current shape in any dimension.  Raises an error if the requested new
    #' shape exceeds `maxshape` in any dimension. Raises an error if the array
    #' doesn't already have a shape: in that case please call
    #' `tiledbsoma_upgrade_shape`.
    #' @param new_shape A vector of integerish, of the same length as the array's `ndim`.
    #' @return No return value
    resize = function(new_shape) {
      stopifnot(
        "resize is not supported for dense arrays until tiledbsoma 1.15" =
          .dense_arrays_can_have_current_domain() || private$.is_sparse,
        "'new_shape' must be a vector of integerish values, of the same length as maxshape" =
          rlang::is_integerish(new_shape, n = self$ndim()) ||
            (bit64::is.integer64(new_shape) && length(new_shape) == self$ndim())
      )
      # Checking slotwise new shape >= old shape, and <= max_shape, is already done in libtiledbsoma
      resize(self$uri, new_shape, .name_of_function(), private$.soma_context)
    },

    #' @description Allows the array to have a resizeable shape as described in the
    #' TileDB-SOMA 1.15 release notes.  Raises an error if the shape exceeds maxshape in any
    #' dimension. Raises an error if the array already has a shape. The methods
    #' `tiledbsoma_upgrade_shape` and `resize` are very similar: the former must be
    #' call on a pre-1.15 array the first time a shape is set on it; the latter must
    #' be used for subsequent resizes on any array which already has upgraded shape.
    #' @param shape A vector of integerish, of the same length as the array's `ndim`.
    #' @return No return value
    tiledbsoma_upgrade_shape = function(shape) {
      stopifnot(
        "tiledbsoma_upgrade_shape is not supported for dense arrays until tiledbsoma 1.15" =
          .dense_arrays_can_have_current_domain() || private$.is_sparse,
        "'shape' must be a vector of integerish values, of the same length as maxshape" =
          rlang::is_integerish(shape, n = self$ndim()) ||
            (bit64::is.integer64(shape) && length(shape) == self$ndim())
      )
      # Checking slotwise new shape >= old shape, and <= max_shape, is already done in libtiledbsoma
      tiledbsoma_upgrade_shape(self$uri, shape, .name_of_function(), private$.soma_context)
    }
  ),
  private = list(
    .is_sparse = NULL,
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
    #  format acceptable for sr_setup and soma_array_reader
    .convert_coords = function(coords) {
      # Ensure coords is a named list, use to select dim points
      stopifnot(
        "'coords' must be a list" = is.list(coords) && length(coords),
        "'coords' must be a list integerish vectors" =
          all(vapply(
            X = coords,
            FUN = function(x) {
              if (is.null(x)) {
                return(TRUE)
              }
              return(
                (is.null(dim(x)) && !is.factor(x)) &&
                  (rlang::is_integerish(x, finite = TRUE) || (bit64::is.integer64(x) && all(is.finite(x)))) &&
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
        "'new_shape' must be an integerish vector with the same length as the array's maxshape" = rlang::is_integerish(new_shape, n = self$ndim(), finite = TRUE) ||
          (bit64::is.integer64(new_shape) && length(new_shape) == self$ndim() && all(is.finite(new_shape)))
      )

      # convert integer to integer64 to match dimension type
      return(bit64::as.integer64(new_shape))
    },

    # Internal marking of one or zero based matrices for iterated reads
    zero_based = NA
  )
)
