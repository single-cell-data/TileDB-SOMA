#' SOMA Collection
#'
#' @description Contains a key-value mapping where the keys are string names and
#' the values are any SOMA-defined foundational or composed type, including
#' \code{\link{SOMACollection}}, \code{\link{SOMADataFrame}},
#' \code{\link{SOMADenseNDArray}}, \code{\link{SOMASparseNDArray}}, or
#' \code{\link{SOMAExperiment}} (lifecycle: maturing).
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
  inherit = SOMACollectionBase
)
