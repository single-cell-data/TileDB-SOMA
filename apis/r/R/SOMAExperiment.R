#' SOMA Experiment
#'
#' @description `SOMAExperiment` is a specialized [`SOMACollection`],
#' representing one or more modes of measurement across a single collection of
#' cells (aka a "multimodal dataset") with fields:
#'
#' - `obs` ([`SOMADataFrame`]): Primary annotations on the observation axis. The
#'   contents of the `soma_joinid` column define the observation index domain
#'   (AKA `obs_id`). All observations for the `SOMAExperiment` must be defined
#'   in this dataframe.
#' - `ms` ([`SOMACollection`] of [`SOMAMeasurement`]): A collection of named
#'   measurements.
#'
#' @export
SOMAExperiment <- R6::R6Class(
  classname = "SOMAExperiment",
  inherit = SOMACollectionBase,

  active = list(
    #' @field Primary annotations on the observation axis. The contents of the
    #' `soma_joinid` column define the observation index domain (AKA `obs_id`).
    #' All observations for the Experiment must be defined in this dataframe.
    obs = function(value) {
      private$get_or_set_soma_field(value, "obs", "SOMADataFrame")
    },

    #' @field A collection of named measurements.
    ms = function(value) {
      private$get_or_set_soma_field(value, "ms", "SOMAMeasurement")
    }
  )
)
