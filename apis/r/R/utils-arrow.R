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
    stop("Unsupported data type: ", x$name, call. = FALSE)
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

#' Create an Arrow field from a TileDB attribute
#' @noRd
arrow_field_from_tiledb_attr <- function(x) {
  stopifnot(inherits(x, "tiledb_attr"))
  arrow::field(
    name = tiledb::name(x),
    type = arrow_type_from_tiledb_type(tiledb::datatype(x)),
    nullable = tiledb::tiledb_attribute_get_nullable(x)
  )
}

#' Create an Arrow schema from a TileDB array schema
#' @noRd
arrow_schema_from_tiledb_schema <- function(x) {
  stopifnot(inherits(x, "tiledb_array_schema"))
  fields <- c(
    lapply(tiledb::dimensions(x), arrow_field_from_tiledb_dim),
    lapply(tiledb::attrs(x), arrow_field_from_tiledb_attr)
  )
  arrow::schema(fields)
}

#' Validate external pointer to ArrowArray
#' @noRd
check_arrow_pointers <- function(arrlst) {
    stopifnot("First argument must be an external pointer to ArrowArray" = check_arrow_array_tag(arrlst[[1]]),
              "Second argument must be an external pointer to ArrowSchema" = check_arrow_schema_tag(arrlst[[2]]))
}
