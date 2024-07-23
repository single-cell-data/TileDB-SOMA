#' SOMA Experiment
#'
#' @description `SOMAExperiment` is a specialized [`SOMACollection`],
#' representing one or more modes of measurement across a single collection of
#' cells (aka a "multimodal dataset") with pre-defined fields: `obs` and `ms`
#' (see _Active Bindings_ below for details). (lifecycle: maturing)
#'
#' @templateVar class SOMAExperiment
#' @template section-add-object-to-collection
 #' @param row_index_name An optional scalar character. If provided, and if
#' the `values` argument is a `data.frame` with row names, then the row
#' names will be extracted and added as a new column to the `data.frame`
#' prior to performing the update. The name of this new column will be set
#' to the value specified by `row_index_name`.
#'
#' @export
SOMAExperiment <- R6::R6Class(
  classname = "SOMAExperiment",
  inherit = SOMACollectionBase,

  public = list(
    #' @description Subset and extract data from a [`SOMAMeasurement`] by
    #' querying the `obs`/`var` axes.
    #' @param measurement_name The name of the measurement to query.
    #' @param obs_query,var_query An [`SOMAAxisQuery`] object for the obs/var
    #' axis.
    #' @return A [`SOMAExperimentAxisQuery`] object.
    axis_query = function(measurement_name, obs_query = NULL, var_query = NULL) {
      SOMAExperimentAxisQuery$new(
        experiment = self,
        measurement_name = measurement_name,
        obs_query = obs_query,
        var_query = var_query
      )
    },

    #' @description Update the obs [`SOMADataFrame`] to add or remove columns.
    #' See [`SOMADataFrame$update()`][`SOMADataFrame`] for details.
    #' @param values A `data.frame`, [`arrow::Table`], or
    #' [`arrow::RecordBatch`].
    update_obs = function(values, row_index_name = NULL) {
      self$obs$update(values, row_index_name)
    },

    #' @description Update the var `SOMADataFrame` to add or remove columns.
    #' See [`SOMADataFrame$update()`][`SOMADataFrame`] for details.
    #' @param values A `data.frame`, [`arrow::Table`], or
    #' [`arrow::RecordBatch`].
    #' @param measurement_name The name of the [`SOMAMeasurement`] whose `var`
    #' will be updated.
    update_var = function(values, measurement_name, row_index_name = NULL) {
      stopifnot(
        "Must specify a single measurement name" =
          is_scalar_character(measurement_name),
        "Measurement does not exist in the experiment" =
          measurement_name %in% self$ms$names()
      )
      self$ms$get(measurement_name)$var$update(values, row_index_name)
    }
  ),

  active = list(
    #' @field obs a [`SOMADataFrame`] containing primary annotations on the
    #' observation axis. The contents of the `soma_joinid` column define the
    #' observation index domain, `obs_id`. All observations for the
    #' `SOMAExperiment` must be defined in this dataframe.
    obs = function(value) {
      private$get_or_set_soma_field(value, "obs", "SOMADataFrame")
    },

    #' @field ms a [`SOMACollection`] of named [`SOMAMeasurement`]s.
    ms = function(value) {
      private$get_or_set_soma_field(value, "ms", "SOMACollection")
    }
  )
)
