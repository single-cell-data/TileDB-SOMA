#' SOMA Axis Query
#' @description Construct a single-axis query object with a combination of
#' coordinates and/or value filters for use with [`SOMAExperimentAxisQuery`].
#' (lifecycle: maturing)
#'
#' Per dimension, the `SOMAAxisQuery` can have value of:
#'
#' * None (i.e., `coords = NULL` and `value_filter = NULL`) - read all values
#' * Coordinates - a set of coordinates on the axis dataframe index, expressed
#'   in any type or format supported by [`SOMADataFrame`]'s `read()` method.
#' * A SOMA `value_filter` across columns in the axis dataframe, expressed as
#'   string
#' * Or, a combination of coordinates and value filter.
#'
#' @seealso [`tiledb::parse_query_condition()`] for more information about valid
#' value filters.
#'
#' @export
SOMAAxisQuery <- R6Class(
  classname = "SOMAAxisQuery",
  public = list(
    #' @field coords The coordinates for the query.
    coords = NULL,
    #' @field value_filter The value filter for the query.
    value_filter = NULL,

    #' @description Create a new `SOMAAxisQuery` object.
    #' @param coords Optional indices specifying the rows to read: either a
    #' vector of the appropriate type or a named list of vectors corresponding
    #' to each dimension.
    #' @param value_filter Optional string containing a logical expression that
    #' is used to filter the returned values.
    initialize = function(value_filter = NULL, coords = NULL) {
      self$coords <- validate_read_coords(coords, dimnames = NULL)
      self$value_filter <- validate_read_value_filter(value_filter)
    }
  )
)
