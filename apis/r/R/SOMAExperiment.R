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
