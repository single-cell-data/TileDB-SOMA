#' SOMA Measurement
#'
#' @description A \code{SOMAMeasurement} is a sub-element of a
#' \code{\link{SOMAExperiment}}, and is otherwise a specialized
#' \code{\link{SOMACollection}} with pre-defined fields: \code{X}, \code{var},
#' \code{obsm}/\code{varm}, and \code{obsp}/\code{varp} (see
#' \emph{Active Bindings} below for details) (lifecycle: maturing).
#'
#' @templateVar class SOMAMeasurement
#' @template section-add-object-to-collection
#'
#' @export
#'
#' @inherit SOMAMeasurementCreate examples
#'
SOMAMeasurement <- R6::R6Class(
  classname = "SOMAMeasurement",
  inherit = SOMACollectionBase,
  active = list(
    #' @field var A \code{\link{SOMADataFrame}} containing primary annotations
    #' on the variable axis, for variables in this measurement
    #' (i.e., annotates columns of \code{X}). The contents of the
    #' \code{soma_joinid} column define the variable index domain,
    #' \code{var_id}. All variables for this measurement must be defined in this
    #' data frame.
    #'
    var = function(value) {
      private$get_or_set_soma_field(value, "var", "SOMADataFrame")
    },

    #' @field X A \code{\link{SOMACollection}} of
    #' \code{\link{SOMASparseNDArray}s}, each contains measured feature values
    #' indexed by \code{[obsid, varid]}.
    #'
    X = function(value) {
      private$get_or_set_soma_field(value, "X", "SOMACollection")
    },

    #' @field obsm A \code{\link{SOMACollection}} of
    #' \code{\link{SOMADenseNDArray}s} containing annotations on the observation
    #' axis. Each array is indexed by \code{obsid} and has the same shape as
    #' \code{obs}.
    #'
    obsm = function(value) {
      private$get_or_set_soma_field(value, "obsm", "SOMACollection")
    },

    #' @field obsp A \code{\link{SOMACollection}} of
    #' \code{\link{SOMASparseNDArray}s} containing pairwise annotations on the
    #' observation axis and indexed with \code{[obsid_1, obsid_2]}.
    #'
    obsp = function(value) {
      private$get_or_set_soma_field(value, "obsp", "SOMACollection")
    },

    #' @field varm A \code{\link{SOMACollection}} of
    #' \code{\link{SOMADenseNDArray}s} containing annotations on the variable
    #' axis. Each array is indexed by \code{varid} and has the same shape as
    #' \code{var}.
    #'
    varm = function(value) {
      private$get_or_set_soma_field(value, "varm", "SOMACollection")
    },

    #' @field varp A \code{\link{SOMACollection}} of
    #' \code{\link{SOMASparseNDArray}s} containing pairwise annotations on the
    #' variable axis and indexed with \code{[varid_1, varid_2]}.
    #'
    varp = function(value) {
      private$get_or_set_soma_field(value, "varp", "SOMACollection")
    }
  )
)
