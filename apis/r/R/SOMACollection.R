#' SOMA Collection
#'
#' @description Contains a key-value mapping where the keys are string names
#' and the values are any SOMA-defined foundational or composed type, including
#' `SOMACollection`, `SOMADataFrame`, `SOMADenseNDArray`, `SOMASparseNDArray`
#' or `SOMAExperiment`.  (lifecycle: experimental)
#'
#' @templateVar id SOMACollection
#' @template section-add-element-to-collection
#'
#' @export
SOMACollection <- R6::R6Class(
  classname = "SOMACollection",
  inherit = SOMACollectionBase
)
