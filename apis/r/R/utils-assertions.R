#' Check if object is empty
#' @noRd
is_empty <- function(x) {
  switch(class(x)[1],
    "data.frame" = nrow(x) == 0,
    length(x) == 0
  )
}

#' Check if a vector is named
#' @noRd
is_named <- function(x, allow_empty = TRUE) {
  !is.null(names(x)) && ifelse(allow_empty, TRUE, all(nzchar(x = names(x = x))))
}

is_named_list <- function(x) {
  is.list(x) && is_named(x)
}

is_scalar_logical <- function(x) {
  is.logical(x) && length(x) == 1
}

is_scalar_character <- function(x) {
  is.character(x) && length(x) == 1
}

is_character_or_null <- function(x) {
  is.character(x) || is.null(x)
}

has_character_rownames <- function(x) {
  stopifnot(is.data.frame(x))
  typeof(attr(x, "row.names")) == "character"
}

is_matrix <- function(x) {
  is.matrix(x) || inherits(x, "Matrix")
}

is_vector_or_int64 <- function(x) {
    is.vector(x) || inherits(x, "integer64")
}

has_dimnames <- function(x) {
  stopifnot(is_matrix(x))
  dims <- dimnames(x) %||% list(NULL)
  all(!vapply(dims, is.null, logical(1L)))
}

string_starts_with <- function(x, prefix) {
  prefix <- paste0("^", prefix)
  grepl(prefix, x)
}

check_package <- function(package) {
  if (requireNamespace(package, quietly = TRUE)) {
    return(invisible())
  }
  stop(paste0("Package '", package, "' must be installed"))
}

#' Assert all values of `x` are a subset of `y`. @param x,y vectors of values
#' @param type A character vector of length 1 used in the error message
#' @return `TRUE` if all values of `x` are present in `y`, otherwise an
#' informative error is thrown with the missing values.
#' @noRd
assert_subset <- function(x, y, type = "value") {
  stopifnot(is.atomic(x) && is.atomic(y))
  missing <- !x %in% y
  if (any(missing)) {
    stop(sprintf(
      "The following %s%s not exist: %s",
      type,
      ifelse(length(missing) == 1, " does", "s do"),
      glue::glue_collapse(x[missing], sep = ", ", last = " and ")
    ), call. = FALSE)
  }
  TRUE
}

#' Validate read coordinates
#'
#' Ensures that coords is one of the following supported formats:
#'
#' - `NULL`
#' - a vector of coordinates to index a single dimension
#' - an unnamed  list of coordinates to index a single dimension
#' - a named list of coordinates to index multiple dimensions (names are
#'   optionally validated against a vector of dimension names)
#'
#' @param coords A vector or list of coordinates
#' @param dimnames character vector of array dimension names
#' @param schema arrow schema object
#' @noRd
validate_read_coords <- function(coords, dimnames = NULL, schema = NULL) {
  # NULL is a valid value
  if (is.null(coords)) return(coords)

  # If coords is a vector, wrap it in a list
  if (is.atomic(coords)) coords <- list(coords)

  # List of multiple coordinate vectors must be named
  if (length(coords) > 1) {
    stopifnot(
      "'coords' must be a named list to query multiple dimensions" =
        is_named_list(coords)
    )
  }

  # TODO: vector type should be validated by array dimension type
  stopifnot(
    "'coords' must be a list of numeric vectors" =
      all(vapply_lgl(coords, is.numeric))
  )

  if (is.null(dimnames)) {
    # Schema isn't useful without dimnames because we don't know which fields
    # are attributes and which are dimensions.
    if (!is.null(schema)) {
      stop(
      "'dimnames' must be provided with a 'schema'", call. = FALSE)
    }
  } else {
    #
    stopifnot(
      "'dimnames' must be a character vector" = is.character(dimnames),
      "names of 'coords' must correspond to dimension names" =
        all(names(coords) %in% dimnames)
    )

    # Apply dimension name to unnamed list containing single coordinate vector
    if (!is_named_list(coords) && length(coords) == 1) {
      coords <- stats::setNames(coords, dimnames)
    }

    # With a schema, we can check if any dimensions are int64 and cast them
    if (!is.null(schema)) {
      stopifnot(
        is_arrow_schema(schema),
        all(dimnames %in% schema$names)
      )

      # identify int64 dimensions
      int64_dims <- vapply_lgl(dimnames, function(x) {
        schema$GetFieldByName(x)$type$ToString() == "int64"
      })

      if (any(int64_dims)) {
        int64_dims <- names(Filter(isTRUE, int64_dims))
        coords[int64_dims] <- recursively_make_integer64(coords[int64_dims])
      }
    }

  }

  coords
}

#' Validate read/query value filter
#' @noRd
validate_read_value_filter <- function(value_filter) {
  stopifnot(
    "'value_filter' must be a scalar character" =
      is.null(value_filter) || is_scalar_character(value_filter)
    )
  value_filter
}

#' Recursively traverse data structure and upcast integer vectors to integer64
#'
#' This is needed as we may receive (named or unnamed) list and/or plain vectors
#' @noRd
recursively_make_integer64 <- function(x) {
    if (is.null(x) || is.character(x) || is.double(x) || is.factor(x) || is.ordered(x)) {
        x 	# do nothing
    } else if (is.list(x)) {
        for (i in seq_along(x)) {
            x[[i]] <- recursively_make_integer64(x[[i]])
        }
    } else if (is.integer(x)) {
        x <- bit64::as.integer64(x)
    } else {
        warning("encountered ", class(x))
    }
    x
}

#' Warn if using a SOMADenseNDArray
#' @param collection_name name of the SOMA collection
#' @param layer a SOMASparseNDArray or SOMADenseNDArray
#' @return Invisible logical indicating if the layer is sparse
#' @references https://github.com/single-cell-data/TileDB-SOMA/issues/1245
#' @noRd

warn_if_dense <- function(collection_name, layer) {
  stopifnot(
    is_scalar_character(collection_name),
    inherits(layer, "SOMASparseNDArray") || inherits(layer, "SOMADenseNDArray")
  )
  is_sparse <- inherits(layer, "SOMASparseNDArray")
  if (!is_sparse) {
    msg1 <- sprintf(
      "Values for '%s' are encoded as dense instead of sparse\n",
      basename(layer$uri)
    )
    msg2 <- sprintf(
      "  - all '%s' arrays should be saved as 'SOMASparseNDArrays'\n",
      collection_name
    )
    warning(msg1, msg2, call. = FALSE, immediate. = TRUE)
  }
  invisible(is_sparse)
}
