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

#' Convert Arrow types to supported TileDB type
#' List of TileDB types supported in R: https://github.com/TileDB-Inc/TileDB-R/blob/8014da156b5fee5b4cc221d57b4aa7d388abc968/inst/tinytest/test_dim.R#L97-L121
#'
#' List of all arrow types: https://github.com/apache/arrow/blob/90aac16761b7dbf5fe931bc8837cad5116939270/r/R/type.R#L700
#' @noRd

tiledb_type_from_arrow_type <- function(x) {
  stopifnot(is_arrow_data_type(x))
  switch(x$name,

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
    # tiledb::r_to_tiledb_type() returns UTF8 for characters but they are
    # not yet queryable so we use ASCII for now
    utf8 = "ASCII",
    string = "ASCII",
    large_utf8 = "ASCII",
    # date32 = "date32",
    # date64 = "date64",
    # time32 = "time32",
    # time64 = "time64",
    # null = "null",
    # timestamp = "timestamp",
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
    dictionary = "INT32", 			# for a dictionary the 'values' are ints, levels are character
    stop("Unsupported Arrow data type: ", x$name, call. = FALSE)
  )
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

#' Retrieve limits for Arrow types
#' @importFrom bit64 lim.integer64
#' @noRd
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
    double =  c(.Machine$double.xmin, .Machine$double.xmax),
    # boolean/bool
    bool = NULL,
    # string/utf8
    utf8 = NULL,
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

#' Create an Arrow field from a TileDB dimension
#' @noRd
arrow_field_from_tiledb_dim <- function(x) {
  stopifnot(inherits(x, "tiledb_dim"))
  arrow::field(
    name = tiledb::name(x),
    type = arrow_type_from_tiledb_type(tiledb::datatype(x)),
    nullable = FALSE
  )
}

## With a nod to Kevin Ushey
#' @noRd
yoink <- function(package, symbol) {
    do.call(":::", list(package, symbol))
}


#' Create an Arrow field from a TileDB attribute
#' @noRd
arrow_field_from_tiledb_attr <- function(x, arrptr=NULL) {
    stopifnot(inherits(x, "tiledb_attr"))
    if (tiledb::tiledb_attribute_has_enumeration(x) && !is.null(arrptr)) {
        .tiledb_array_is_open <- yoink("tiledb", "libtiledb_array_is_open")
        if (!.tiledb_array_is_open(arrptr)) {
            .tiledb_array_open_with_ptr <- yoink("tiledb", "libtiledb_array_open_with_ptr")
            arrptr <- .tiledb_array_open_with_ptr(arrptr, "READ")
        }
        ord <- tiledb::tiledb_attribute_is_ordered_enumeration_ptr(x, arrptr)
        idx <- arrow_type_from_tiledb_type(tiledb::datatype(x))
        arrow::field(name = tiledb::name(x),
                     type = arrow::dictionary(index_type=idx, ordered=ord),
                     nullable = tiledb::tiledb_attribute_get_nullable(x))
    } else {
        arrow::field(name = tiledb::name(x),
                     type = arrow_type_from_tiledb_type(tiledb::datatype(x)),
                     nullable = tiledb::tiledb_attribute_get_nullable(x))
    }
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

  field_type <- tiledb_type_from_arrow_type(field$type)
  tiledb::tiledb_attr(
    name = field$name,
    type = field_type,
    nullable = field$nullable,
    ncells = if (field_type == "ASCII") NA_integer_ else 1L,
    filter_list = tiledb::tiledb_filter_list(
      tiledb_create_options$attr_filters(
        attr_name = field$name,
        default = list(default_zstd_filter)
      )
    )
  )
}

#' Create an Arrow schema from a TileDB array schema
#' @noRd
arrow_schema_from_tiledb_schema <- function(x) {
  stopifnot(inherits(x, "tiledb_array_schema"))
  dimfields <- lapply(tiledb::dimensions(x), arrow_field_from_tiledb_dim)
  if (!is.null(x@arrptr)) {
      attfields <- lapply(tiledb::attrs(x), arrow_field_from_tiledb_attr, x@arrptr)
  } else {
      attfields <- lapply(tiledb::attrs(x), arrow_field_from_tiledb_attr)
  }
  arrow::schema(c(dimfields, attfields))
}

#' Validate external pointer to ArrowArray
#' @noRd
check_arrow_pointers <- function(arrlst) {
    stopifnot("First argument must be an external pointer to ArrowArray" = check_arrow_array_tag(arrlst[[1]]),
              "Second argument must be an external pointer to ArrowSchema" = check_arrow_schema_tag(arrlst[[2]]))
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
    "'from' and 'to' must both be Arrow DataTypes"
      = is_arrow_data_type(from) && is_arrow_data_type(to)
  )

  is_string <- function(x) {
    x$ToString() %in% c("string", "large_string")
  }

  compatible <- if (is_string(from) && is_string(to)) {
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
    "'from' and 'to' must both be Arrow Schemas"
      = is_arrow_schema(from) && is_arrow_schema(to),
    "'from' and 'to' must have the same number of fields"
      = length(from) == length(to),
    "'from' and 'to' must have the same field names"
      = identical(sort(names(from)), sort(names(to)))
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
extract_levels <- function(arrtbl, exclude_cols=c("soma_joinid")) {
    stopifnot("Argument must be an Arrow Table object" = is_arrow_table(arrtbl))
    nm <- names(arrtbl)                 # we go over the table column by column
    nm <- nm[-match(exclude_cols, nm)]  # but skip soma_joinid etc as in exclude_cols
    reslst <- vector(mode = "list", length = length(nm))
    names(reslst) <- nm		# and fill a named list, entries default to NULL
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
