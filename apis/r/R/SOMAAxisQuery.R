#' @description Single-axis dataframe query with coordinates and a value filter.
#' [lifecycle: experimental]
#' Per dimension, the AxisQuery can have value of:
#'
#' * None - all data
#' * Coordinates - a set of coordinates on the axis dataframe index, expressed
#'   in any type or format supported by ``DataFrame.read()``.
#' * A SOMA ``value_filter`` across columns in the axis dataframe, expressed as
#'   string
#' * Or, a combination of coordinates and value filter.
#'
#' @export
AxisQuery <- R6Class(
  classname = "AxisQuery",
  public = list(
    value_filter = NULL,
    coords = NULL,

    initialize = function(value_filter = NULL, coords = NULL) {
      self$coords <- validate_read_coords(coords, dimnames = NULL)
      self$value_filter <- validate_read_value_filter(value_filter)
    }
  )
)

#' Validate read coordinates
#' Ensures that coords is a named list of numeric vectors and, optionally, that
#' the names of the list correspond to the dimension names of the array.
#' @param coords A list of coordinates
#' @param dimnames vector of array dimension names
#' @noRd
validate_read_coords <- function(coords, dimnames = NULL) {
  stopifnot(
    "'coords' must be a named list" =
      is.null(coords) || is_named_list(coords),
    "'coords' must be a list of numeric vectors" =
      all(vapply_lgl(coords, is.numeric))
  )

  if (!is.null(dimnames)) {
    stopifnot(
      "'dimnames' must be a character vector" = is.character(dimnames),
      "names of 'coords' must correspond to dimension names" =
        all(names(coords) %in% dimnames)
    )
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
