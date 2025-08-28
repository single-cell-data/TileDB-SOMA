is_arrow_object <- function(x) {
  inherits(x, "ArrowObject")
}

is_arrow_data_type <- function(x) {
  is_arrow_object(x) && inherits(x, "DataType")
}

is_arrow_field <- function(x) {
  is_arrow_object(x) && inherits(x, "Field")
}

is_arrow_record_batch <- function(x) {
  is_arrow_object(x) && inherits(x, "RecordBatch")
}

is_arrow_array <- function(x) {
  is_arrow_object(x) && inherits(x, "Array")
}

is_arrow_chunked_array <- function(x) {
  is_arrow_object(x) && inherits(x, "ChunkedArray")
}

is_arrow_table <- function(x) {
  is_arrow_object(x) && inherits(x, "Table")
}

is_arrow_schema <- function(x) {
  is_arrow_object(x) && inherits(x, "Schema")
}

is_arrow_dictionary <- function(x) {
  is_arrow_object(x) && inherits(x, "Field") && inherits(x$type, "DictionaryType")
}

#' @method as.logical Scalar
#' @export
#'
as.logical.Scalar <- \(x, ...) as.logical(x$as_vector(), ...)

#' Convert Arrow types to supported TileDB type
#' List of TileDB types supported in R: https://github.com/TileDB-Inc/TileDB-R/blob/8014da156b5fee5b4cc221d57b4aa7d388abc968/inst/tinytest/test_dim.R#L97-L121
#' Note: TileDB attrs may be UTF-8; TileDB dims may not.
#'
#' List of all arrow types: https://github.com/apache/arrow/blob/90aac16761b7dbf5fe931bc8837cad5116939270/r/R/type.R#L700
#' @noRd

tiledb_type_from_arrow_type <- function(x, is_dim) {
  stopifnot(is_arrow_data_type(x))
  retval <- switch(x$name,
    int8 = "INT8",
    int16 = "INT16",
    int32 = "INT32",
    int64 = "INT64",
    uint8 = "UINT8",
    uint16 = "UINT16",
    uint32 = "UINT32",
    uint64 = "UINT64",
    float32 = "FLOAT32",
    float = "FLOAT32",
    float64 = "FLOAT64",
    # based on tiledb::r_to_tiledb_type()
    double = "FLOAT64",
    boolean = "BOOL",
    bool = "BOOL",
    # large_string = "large_string",
    # binary = "binary",
    # large_binary = "large_binary",
    # fixed_size_binary = "fixed_size_binary",
    utf8 = "UTF8",
    string = "UTF8",
    large_utf8 = "UTF8",
    # based on what TileDB supports
    date32 = "DATETIME_DAY",
    # date64 = "date64",
    # time32 = "time32",
    # time64 = "time64",
    # null = "null",
    # based on what TileDB supports with a default msec res.
    timestamp = "DATETIME_MS",
    # decimal128 = "decimal128",
    # decimal256 = "decimal256",
    # struct = "struct",
    # list_of = "list",
    # list = "list",
    # large_list_of = "large_list",
    # large_list = "large_list",
    # fixed_size_list_of = "fixed_size_list",
    # fixed_size_list = "fixed_size_list",
    # map_of = "map",
    # duration = "duration",
    dictionary = tiledb_type_from_arrow_type(x$index_type, is_dim = is_dim),
    stop("Unsupported Arrow data type: ", x$name, call. = FALSE)
  )
  if (is_dim && retval == "UTF8") {
    retval <- "ASCII"
  }
  retval
}

arrow_type_from_tiledb_type <- function(x) {
  stopifnot(is.character(x))
  switch(x,
    INT8 = arrow::int8(),
    INT16 = arrow::int16(),
    INT32 = arrow::int32(),
    INT64 = arrow::int64(),
    UINT8 = arrow::uint8(),
    UINT16 = arrow::uint16(),
    UINT32 = arrow::uint32(),
    UINT64 = arrow::uint64(),
    FLOAT32 = arrow::float32(),
    FLOAT64 = arrow::float64(),
    BOOL = arrow::boolean(),
    ASCII = arrow::utf8(),
    UTF8 = arrow::utf8(),
    stop("Unsupported data type: ", x, call. = FALSE)
  )
}

#' Get the \R Type from an Arrow Type
#'
#' Get an \R \link[base:typeof]{type} from an Arrow type. This function is
#' equivalent to \code{\link[base]{typeof}()} rather than
#' \code{\link[base]{mode}()} or \code{\link[base]{class}()}, and returns the
#' equivalent \strong{type}. For example, the equivalent \R type to an Arrow
#' \link[arrow]{dictionary} is \dQuote{\code{integer}}, not
#' \dQuote{\code{factor}}; likewise, the equivalent \R type to an Arrow 64-bit
#' integer is \dQuote{\code{double}}
#'
#' @param x An \CRANpkg{arrow} \link[arrow:Schema]{schema},
#' \link[arrow:Field]{field}, or \link[arrow:infer_type]{data type}
#'
#' @return If \code{x} is a \link[arrow:infer_type]{data type}, a single
#' character value giving the \R \link[base:typeof]{type} of \code{x}; if no
#' corresponding \R type, returns the \CRANpkg{arrow} type name
#'
#' @return If \code{x} is a \link[arrow:Field]{field}, a single named character
#' vector with the name being the field name and the value being the \R
#' \link[base:typeof]{type}
#'
#' @return If \code{x} is a \link[arrow:Schema]{schema}, a named vector where
#' the names are field names and the values are the \R \link[base:typeof]{types}
#' of each field
#'
#' @keywords internal
#'
# @export
#'
#' @seealso \code{\link[base]{typeof}()}
#'
#' @noRd
#'
r_type_from_arrow_type <- function(x) UseMethod("r_type_from_arrow_type")

#' @rdname r_type_from_arrow_type
#'
#' @method r_type_from_arrow_type Schema
#' @export
#'
#' @noRd
#'
r_type_from_arrow_type.Schema <- function(x) {
  return(vapply(
    X = x$names,
    FUN = function(f) r_type_from_arrow_type(x[[f]]),
    FUN.VALUE = character(1L),
    USE.NAMES = TRUE
  ))
}

#' @rdname r_type_from_arrow_type
#'
#' @method r_type_from_arrow_type Field
#' @export
#'
#' @noRd
#'
r_type_from_arrow_type.Field <- function(x) {
  tt <- r_type_from_arrow_type(x$type)
  names(x = tt) <- x$name
  return(tt)
}

#' @rdname r_type_from_arrow_type
#'
#' @method r_type_from_arrow_type DataType
#' @export
#'
#' @noRd
#'
r_type_from_arrow_type.DataType <- function(x) {
  # Types are equivalent to `typeof()`, not `mode()` or `class()`
  return(switch(
    EXPR = x$name,
    int8 = ,
    int16 = ,
    int32 = ,
    dictionary = ,
    uint8 = ,
    uint16 = ,
    uint32 = "integer",
    int64 = ,
    uint64 = ,
    date32 = ,
    timestamp = ,
    float = "double",
    bool = "logical",
    utf8 = ,
    large_utf8 = "character",
    x$name
  ))
}

#' Retrieve limits for Arrow types
#'
#' @noRd
#'
arrow_type_range <- function(x) {
  stopifnot(is_arrow_data_type(x))

  switch(x$name,
    int8 = c(-128L, 127L),
    int16 = c(-32768L, 32767L),
    int32 = c(-2147483647L, 2147483647L),
    int64 = bit64::lim.integer64(),
    uint8 = c(0L, 255L),
    uint16 = c(0L, 65535L),
    uint32 = bit64::as.integer64(c(0, 4294967295)),
    # We can't specify the full range of uint64 in R so we use the max of int64
    uint64 = c(bit64::as.integer64(0), bit64::lim.integer64()[2]),
    # float32/float
    float = c(-3.4028235e+38, 3.4028235e+38),
    # float64/double
    double = c(.Machine$double.xmin, .Machine$double.xmax),
    # boolean/bool
    bool = NULL,
    # string/utf8
    utf8 = NULL,
    large_utf8 = NULL,
    stop("Unsupported data type:", x$name, call. = FALSE)
  )
}

#' Retrieve unsigned limits for Arrow types
#' This restricts the lower bound of signed numeric types to 0
#' @noRd
arrow_type_unsigned_range <- function(x) {
  range <- arrow_type_range(x)
  range[1] <- switch(x$name,
    int8 = 0L,
    int16 = 0L,
    int32 = 0L,
    int64 = bit64::as.integer64(0),
    float = 0,
    double = 0,
    range[1]
  )
  range
}

#' Create a TileDB attribute from an Arrow field
#' @return a [`tiledb::tiledb_attr-class`]
#' @noRd
tiledb_attr_from_arrow_field <- function(field, tiledb_create_options) {
  stopifnot(
    is_arrow_field(field),
    inherits(tiledb_create_options, "TileDBCreateOptions")
  )

  # Default zstd filter to use if none is specified in platform config
  default_zstd_filter <- list(
    name = "ZSTD",
    COMPRESSION_LEVEL = tiledb_create_options$dataframe_dim_zstd_level()
  )

  field_type <- tiledb_type_from_arrow_type(field$type, is_dim = FALSE)
  tiledb::tiledb_attr(
    name = field$name,
    type = field_type,
    nullable = field$nullable,
    ncells = if (field_type == "ASCII" || field_type == "UTF8") NA_integer_ else 1L,
    filter_list = tiledb::tiledb_filter_list(
      tiledb_create_options$attr_filters(
        attr_name = field$name,
        default = list(default_zstd_filter)
      )
    )
  )
}

#' Validate external pointer to ArrowArray which is embedded in a nanoarrow S3 type
#' @noRd
check_arrow_pointers <- function(arrlst) {
  stopifnot(inherits(arrlst, "nanoarrow_array"))
}

#' Validate compatibility of Arrow data types
#'
#' For most data types, this is a simple equality check but it also provides
#' allowances for certain comparisons:
#'
#' - string and large_string
#'
#' @param from an [`arrow::DataType`]
#' @param to an [`arrow::DataType`]
#' @return a logical indicating whether the data types are compatible
#' @noRd
check_arrow_data_types <- function(from, to) {
  stopifnot(
    "'from' and 'to' must both be Arrow DataTypes" = is_arrow_data_type(from) && is_arrow_data_type(to)
  )

  is_string <- function(x) {
    x$ToString() %in% c("string", "large_string")
  }

  is_dict_of_string <- function(x) {
    startsWith(x$ToString(), "dictionary<values=large_string,") ||
      startsWith(x$ToString(), "dictionary<values=string,")
  }

  compatible <- if ((is_string(from) && is_string(to)) || (is_dict_of_string(from) && is_dict_of_string(to))) {
    TRUE
  } else {
    from$Equals(to)
  }

  compatible
}

#' Validate compatibility of Arrow schemas
#'
#' This is essentially a vectorized version of [`check_arrow_data_types`] that
#' checks the compatibility of each field in the schemas.
#' @param from an [`arrow::Schema`]
#' @param to an [`arrow::Schema`] with the same set of fields as `from`
#' @return `TRUE` if the schemas are compatible, otherwise an error is thrown
#' @noRd
check_arrow_schema_data_types <- function(from, to) {
  stopifnot(
    "'from' and 'to' must both be Arrow Schemas" = is_arrow_schema(from) && is_arrow_schema(to),
    "'from' and 'to' must have the same number of fields" = length(from) == length(to),
    "'from' and 'to' must have the same field names" = identical(sort(names(from)), sort(names(to)))
  )

  fields <- names(from)
  msgs <- character(0L)
  for (field in fields) {
    from_type <- from[[field]]$type
    to_type <- to[[field]]$type
    if (!check_arrow_data_types(from_type, to_type)) {
      msg <- sprintf(
        "  - field '%s': %s != %s\n",
        field,
        from_type$ToString(),
        to_type$ToString()
      )
      msgs <- c(msgs, msg)
    }
  }

  if (length(msgs) > 0L) {
    stop(
      "Schemas are incompatible:\n",
      string_collapse(msgs, sep = "\n"),
      call. = FALSE
    )
  }
  return(TRUE)
}

#' Extract levels from dictionaries
#' @noRd
extract_levels <- function(arrtbl, exclude_cols = c("soma_joinid")) {
  stopifnot("Argument must be an Arrow Table object" = is_arrow_table(arrtbl))
  nm <- names(arrtbl) # we go over the table column by column
  nm <- nm[-match(exclude_cols, nm)] # but skip soma_joinid etc as in exclude_cols
  reslst <- vector(mode = "list", length = length(nm))
  names(reslst) <- nm # and fill a named list, entries default to NULL
  for (n in nm) {
    inftp <- arrow::infer_type(arrtbl[[n]])
    if (inherits(inftp, "DictionaryType")) {
      # levels() extracts the enumeration levels from the factor vector we have
      reslst[[n]] <- levels(arrtbl[[n]]$as_vector())
      # set 'ordered' attribute
      attr(reslst[[n]], "ordered") <- inftp$ordered
    }
  }
  reslst
}


#' Domain and extent table creation helper for data.frame writes returning a Table with
#' a column per dimension for the given (incoming) arrow schema of a Table
#' @noRd
get_domain_and_extent_dataframe <- function(
  tbl_schema,
  ind_col_names,
  domain = NULL,
  tdco = TileDBCreateOptions$new(PlatformConfig$new())
) {
  stopifnot(
    "First argument must be an arrow schema" = inherits(tbl_schema, "Schema"),
    "Second argument must be character" = is.character(ind_col_names),
    "Second argument cannot be empty vector" = length(ind_col_names) > 0,
    "Second argument index names must be columns in first argument" =
      all(is.finite(match(ind_col_names, names(tbl_schema)))),
    "Third argument must be options wrapper" = inherits(tdco, "TileDBCreateOptions")
  )
  stopifnot(
    "domain must be NULL or a named list, with values being 2-element vectors or NULL" = is.null(domain) ||
      ( # Check that `domain` is a list of length `length(ind_col_names)`
        # where all values are named after `ind_col_names`
        # and all values are `NULL` or a two-length atomic non-factor vector
        rlang::is_list(domain, n = length(ind_col_names)) &&
          identical(sort(names(domain)), sort(ind_col_names)) &&
          all(vapply_lgl(
            domain,
            function(x) is.null(x) || (is.atomic(x) && !is.factor(x) && length(x) == 2L)
          ))
      )
  )

  rl <- sapply(ind_col_names, \(ind_col_name) {
    ind_col <- tbl_schema$GetFieldByName(ind_col_name)
    ind_col_type <- ind_col$type
    ind_col_type_name <- ind_col$type$name

    ind_ext <- tdco$dim_tile(ind_col_name)

    # Default 2048 mods to 0 for 8-bit types and 0 is an invalid extent
    if (ind_col$type$bit_width %||% 0L == 8L) {
      ind_ext <- 64L
    }

    # We need to subtract off extent from the max because if we don't:
    #
    # Error: [TileDB::Dimension] Error: Tile extent check failed; domain max
    # expanded to multiple of tile extent exceeds max value representable by
    # domain type. Reduce domain max by 1 tile extent to allow for
    # expansion.
    if (ind_col_name == "soma_joinid") {
      # Must be non-negative
      ind_max_dom <- arrow_type_unsigned_range(ind_col_type) - c(0, ind_ext)
    } else {
      # Others can be negative
      ind_max_dom <- arrow_type_range(ind_col_type) - c(0, ind_ext)
    }

    requested_slot <- domain[[ind_col_name]]
    ind_cur_dom <- if (is.null(requested_slot)) {
      # New shape: if the slot is null, make the size as small
      # as possible since current domain can only be resized upward.
      #
      # Core current-domain semantics are (lo, hi) with both
      # inclusive, with lo <= hi. This means smallest is (0, 0)
      # which is shape 1, not 0.
      if (bit64::is.integer64(ind_max_dom)) {
        c(bit64::as.integer64(0), bit64::as.integer64(0))
      } else if (is.integer(ind_max_dom)) {
        c(0L, 0L)
      } else {
        c(0, 0)
      }
    } else {
      requested_slot
    }
    # Core supports no domain specification for variable-length dims, which
    # includes string/binary dims.
    if (ind_col_type_name %in% c("string", "large_utf8", "utf8")) ind_ext <- NA

    # https://github.com/single-cell-data/TileDB-SOMA/issues/2407
    if (ind_col_type_name %in% c("string", "utf8", "large_utf8")) {
      aa <- if (is.null(requested_slot)) {
        arrow::arrow_array(c("", "", "", "", ""), ind_col_type)
      } else {
        arrow::arrow_array(c("", "", "", requested_slot[[1]], requested_slot[[2]]), ind_col_type)
      }
    } else {
      # If they wanted (0, 99) then extent must be at most 100.
      # This is tricky though. Some cases:
      # * lo = 0, hi = 99, extent = 1000
      #   We look at hi - lo + 1; resize extent down to 100
      # * lo = 1000, hi = 1099, extent = 1000
      #   We look at hi - lo + 1; resize extent down to 100
      # * lo = min for datatype, hi = max for datatype
      #   We get integer overflow trying to compute hi - lo + 1
      # So if lo <= 0 and hi >= ind_ext, this is fine without
      # computing hi - lo + 1.
      lo <- ind_max_dom[[1]]
      hi <- ind_max_dom[[2]]
      if (lo > 0 || hi < ind_ext) {
        dom_span <- hi - lo + 1
        if (ind_ext > dom_span) {
          ind_ext <- dom_span
        }
      }
      aa <- arrow::arrow_array(c(ind_max_dom, ind_ext, ind_cur_dom), ind_col_type)
    }

    aa
  })
  names(rl) <- ind_col_names
  dom_ext_tbl <- do.call(arrow::arrow_table, rl)
  dom_ext_tbl
}

#' Domain and extent table creation helper for array writes returning a Table with
#' a column per dimension for the given (incoming) arrow schema of a Table
#' @noRd
get_domain_and_extent_array <- function(shape, is_sparse) {
  stopifnot("First argument must be vector of positive values" = is.vector(shape) && all(shape > 0))
  indvec <- seq_len(length(shape)) - 1 # sequence 0, ..., length()-1
  rl <- sapply(indvec, \(ind) {
    ind_col <- sprintf("soma_dim_%d", ind)
    ind_col_type <- arrow::int64()

    # TODO:  this function needs to take a
    # TileDBCreateOptions$new(PlatformConfig option as
    # get_domain_and_extent_dataframe does.
    # https://github.com/single-cell-data/TileDB-SOMA/issues/2966
    # For now, the core extent is not taken from the platform_config.
    ind_ext <- shape[ind + 1]

    ind_cur_dom <- c(0L, shape[ind + 1] - 1L)

    # We need to do this because if we don't:
    #
    # Error: [TileDB::Dimension] Error: Tile extent check failed; domain max
    # expanded to multiple of tile extent exceeds max value representable by
    # domain type. Reduce domain max by 1 tile extent to allow for
    # expansion.
    ind_max_dom <- arrow_type_unsigned_range(ind_col_type) - c(0, ind_ext)

    aa <- arrow::arrow_array(c(ind_max_dom, ind_ext, ind_cur_dom), ind_col_type)

    aa
  })
  names(rl) <- sprintf("soma_dim_%d", indvec)
  dom_ext_tbl <- do.call(arrow::arrow_table, rl)
  dom_ext_tbl
}
