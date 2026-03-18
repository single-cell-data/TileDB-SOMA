#' SOMA Collection
#'
#' @description Contains a key-value mapping where the keys are string names and
#' the values are any SOMA-defined foundational or composed type, including
#' \code{\link{SOMACollection}}, \code{\link{SOMADataFrame}},
#' \code{\link{SOMADenseNDArray}}, \code{\link{SOMASparseNDArray}}, or
#' \code{\link{SOMAExperiment}} (lifecycle: maturing).
#'
#' @inherit SOMACollectionBase details
#'
#' @templateVar class SOMACollection
#' @template section-add-object-to-collection
#'
#' @export
#'
#' @inherit SOMACollectionCreate examples
#'
SOMACollection <- R6::R6Class(
  classname = "SOMACollection",
  inherit = SOMACollectionBase,
  private = list(
    # @description Open the handle for the C++ interface
    .open_handle = function(open_mode, timestamp) {
      private$.set_handle(open_collection_handle(
        self$uri,
        open_mode,
        private$.context$handle,
        timestamp
      ))
    },

    # @description Implementation for creating a collection.
    .create = function() {
      soma_collection_create(
        uri = self$uri,
        context = private$.context$handle,
        timestamp = self$.tiledb_timestamp_range
      )
    }
  )
)
