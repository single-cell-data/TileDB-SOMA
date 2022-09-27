is_arrow_object <- function(x) {
  inherits(x, "ArrowObject")
}

is_arrow_data_type <- function(x) {
  is_arrow_object(x) && inherits(x, "DataType")
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
    utf8 = "UTF8",
    # large_utf8 = "large_string",
    # large_string = "large_string",
    # binary = "binary",
    # large_binary = "large_binary",
    # fixed_size_binary = "fixed_size_binary",
    # based on tiledb::r_to_tiledb_type()
    string = "UTF8",
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
    stop("Unsupported data type", call. = FALSE)
  )
}
