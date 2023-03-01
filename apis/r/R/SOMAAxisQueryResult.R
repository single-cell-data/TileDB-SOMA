#' `SOMAExperiment` Axis Query Result
#' @description Access [`SOMAExperimentAxisQuery`] results.
#' @export
SOMAAxisQueryResult <- R6Class(
  classname = "SOMAAxisQueryResult",
  public = list(

    #' @description Create a new `SOMAAxisQueryResult` object.
    #' @param obs,var [`arrow::Table`] containing `obs` or `var` query slice.
    #' @param X_layers named list of [`arrow::Table`]s, one for each `X` layer.
    initialize = function(obs, var, X_layers) {
      stopifnot(
        is_arrow_table(obs),
        is_arrow_table(var),
        is_named_list(X_layers),
        all(vapply_lgl(X_layers, is_arrow_table))
      )

      private$.obs <- obs
      private$.var <- var
      private$.X_layers <- X_layers
    }
  ),

  active = list(
    #' @field obs [`arrow::Table`] containing `obs` query slice.
    obs = function(value) {
      if (!missing(value)) read_only_error("obs")
      private$.obs
    },

    #' @field var [`arrow::Table`] containing `var` query slice.
    #' `measurement_name`.
    var = function(value) {
      if (!missing(value)) read_only_error("var")
      private$.var
    },

    #' @field X_layers named list of [`arrow::Table`]s for each `X` layer.
    X_layers = function(value) {
      if (!missing(value)) read_only_error("ms")
      private$.X_layers
    }
  ),

  private = list(
    .obs = NULL,
    .var = NULL,
    .X_layers = NULL
  )
)
