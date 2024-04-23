#' `SOMAExperiment` Axis Query Result Indexer
#' @description Index `obs`/`var` soma_joinids for a given query result.
#'
#' Retrieve the index of the given `obs` or `var` coordinates in the query
#' result. Coordinates outside of the query result will return
#' [`arrow::null()`].
#' @keywords internal
#' @export
SOMAAxisIndexer <- R6::R6Class("SOMAAxisIndexer",
  public = list(

    #' @description Create a new `SOMAAxisIndexer` object.
    #' @param query The [`SOMAExperimentAxisQuery`] object to build indices for.
    initialize = function(query) {
      stopifnot(inherits(query, "SOMAExperimentAxisQuery"))
      private$.query <- query
    },

    #' @description Get the index of the given `obs` coordinates.
    #' @param coords vector or [`arrow::Array`] of numeric coordinates.
    by_obs = function(coords) {
      if (is.null(coords)) {
        stop("COORDS IS NULL")
      }
      arrow::match_arrow(
        x = private$.validate_coords(coords),
        table = private$.obs_index()
      )
    },

    #' @description Get the index of the given `var` coordinates.
    #' @param coords vector or [`arrow::Array`] of numeric coordinates.
    by_var = function(coords) {
      arrow::match_arrow(
        x = private$.validate_coords(coords),
        table = private$.var_index()
      )
    }
  ),

  private = list(
    .cached_obs = NULL,
    .cached_var = NULL,
    .query = NULL,

    # Retrieve index for the obs axis
    .obs_index = function() {
      if (is.null(private$.cached_obs)) {
        private$.cached_obs <- private$.query$obs_joinids()
      }
      private$.cached_obs
    },

    # Retrieve index for the var axis
    .var_index = function() {
      if (is.null(private$.cached_var)) {
        private$.cached_var <- private$.query$var_joinids()
      }
      private$.cached_var
    },

    .validate_coords = function(coords) {
      print("HEY");
      print("COORDS<<");
      print(coords);
      print(">>COORDS");
      stopifnot(
        "'coords' must be a numeric vector or arrow Array" =
          is.numeric(coords) || is_arrow_array(coords) || is_arrow_chunked_array(coords)
      )
      coords
    }
  )
)
