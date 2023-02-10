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
