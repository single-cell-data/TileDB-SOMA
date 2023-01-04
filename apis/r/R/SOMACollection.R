#' SOMA Collection
#'
#' @description Contains a key-value mapping where the keys are string names
#' and the values are any SOMA-defined foundational or composed type, including
#' `SOMACollection`, `SOMADataFrame`, `SOMADenseNdArray`, `SOMASparseNdArray`
#' or `SOMAExperiment`.
#'
#' @export
SOMACollection <- R6::R6Class(
  classname = "SOMACollection",
  inherit = SOMACollectionBase
)
