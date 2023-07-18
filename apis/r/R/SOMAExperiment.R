#' SOMA Experiment
#'
#' @description `SOMAExperiment` is a specialized [`SOMACollection`],
#' representing one or more modes of measurement across a single collection of
#' cells (aka a "multimodal dataset") with pre-defined fields: `obs` and `ms`
#' (see _Active Bindings_ below for details). (lifecycle: experimental)
#'
#' @templateVar class SOMAExperiment
#' @template section-add-object-to-collection
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
