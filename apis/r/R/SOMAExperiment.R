#' SOMA Experiment
#'
#' @description \code{SOMAExperiment} is a specialized
#' \code{\link{SOMACollection}}, representing one or more modes of measurement
#' across a single collection of cells (aka a \dQuote{multimodal dataset}) with
#' pre-defined fields: \code{obs} and \code{ms} (see \emph{Active Bindings} below
#' for details) (lifecycle: maturing).
#'
#' @templateVar class SOMAExperiment
#' @template section-add-object-to-collection
#'
#' @param values A data frame, \link[arrow:Table]{Arrow table}, or
#' \link[arrow:RecordBatch]{Arrow record batch}.
#' @param row_index_name An optional scalar character. If provided, and if
#' the \code{values} argument is a data frame with row names, then the row
#' names will be extracted and added as a new column to the data frame
#' prior to performing the update. The name of this new column will be set
#' to the value specified by \code{row_index_name}.
#'
#' @export
#'
#' @inherit SOMAExperimentCreate examples
#'
SOMAExperiment <- R6::R6Class(
  classname = "SOMAExperiment",
  inherit = SOMACollectionBase,
  public = list(
    #' @description Subset and extract data from a \code{\link{SOMAMeasurement}}
    #' by querying the \code{obs}/\code{var} axes.
    #'
    #' @param measurement_name The name of the measurement to query.
    #' @param obs_query,var_query An \code{\link{SOMAAxisQuery}} object for the
    #' obs/var axis.
    #'
    #' @return A \code{\link{SOMAExperimentAxisQuery}} object.
    #'
    axis_query = function(measurement_name, obs_query = NULL, var_query = NULL) {
      SOMAExperimentAxisQuery$new(
        experiment = self,
        measurement_name = measurement_name,
        obs_query = obs_query,
        var_query = var_query
      )
    },

    #' @description Update the \code{obs} data frame to add or remove columns.
    #' See \code{\link[tiledbomsa:SOMADataFrame]{SOMADataFrame$update()}} for
    #' more details.
    #'
    update_obs = function(values, row_index_name = NULL) {
      self$obs$update(values, row_index_name)
    },

    #' @description Update the \code{var} data frame to add or remove columns.
    #' See \code{\link[tiledbomsa:SOMADataFrame]{SOMADataFrame$update()}} for
    #' more details.
    #'
    #' @param measurement_name The name of the \code{\link{SOMAMeasurement}}
    #' whose \code{var} will be updated.
    #'
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
    #' @field obs A \code{\link{SOMADataFrame}} containing primary annotations
    #' on the observation axis. The contents of the \code{soma_joinid} column
    #' define the observation index domain, \code{obs_id}. All observations for
    #' the \code{SOMAExperiment} must be defined in this data frame.
    #'
    obs = function(value) {
      private$get_or_set_soma_field(value, "obs", "SOMADataFrame")
    },

    #' @field ms A \code{\link{SOMACollection}} of named
    #' \code{\link{SOMAMeasurement}s}.
    #'
    ms = function(value) {
      private$get_or_set_soma_field(value, "ms", "SOMACollection")
    }
  )
)
