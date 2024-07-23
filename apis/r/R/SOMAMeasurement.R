#' SOMA Measurement
#'
#' @description A `SOMAMeasurement` is a sub-element of a [`SOMAExperiment`],
#' and is otherwise a specialized [`SOMACollection`] with pre-defined fields:
#' `X`, `var`, `obsm`/`varm`, and `obsp`/`varp` (see _Active Bindings_ below for
#' details). (lifecycle: maturing)
#'
#' @templateVar class SOMAMeasurement
#' @template section-add-object-to-collection
#'
#' @export
SOMAMeasurement <- R6::R6Class(
  classname = "SOMAMeasurement",
  inherit = SOMACollectionBase,

  active = list(
    #' @field var a [`SOMADataFrame`] containing primary annotations on the
    #' variable axis, for variables in this measurement (i.e., annotates columns
    #' of `X`). The contents of the `soma_joinid` column define the variable
    #' index domain, `var_id`. All variables for this measurement must be
    #' defined in this dataframe.
    var = function(value) {
      private$get_or_set_soma_field(value, "var", "SOMADataFrame")
    },

    #' @field X a [`SOMACollection`] of [`SOMASparseNDArray`]s, each contains
    #' measured feature values indexed by `[obsid, varid]`.
    X = function(value) {
      private$get_or_set_soma_field(value, "X", "SOMACollection")
    },

    #' @field obsm a [`SOMACollection`] of [`SOMADenseNDArray`]s containing
    #' annotations on the observation axis. Each array is indexed by `obsid` and
    #' has the same shape as `obs`.
    obsm = function(value) {
      private$get_or_set_soma_field(value, "obsm", "SOMACollection")
    },

    #' @field obsp a [`SOMACollection`] of [`SOMASparseNDArray`]s containing
    #' pairwise annotations on the observation axis and indexed with `[obsid_1,
    #' obsid_2]`.
    obsp = function(value) {
      private$get_or_set_soma_field(value, "obsp", "SOMACollection")
    },

    #' @field varm a [`SOMACollection`] of [`SOMADenseNDArray`]s containing
    #' annotations on the variable axis. Each array is indexed by `varid` and
    #' has the same shape as `var`.
    varm = function(value) {
      private$get_or_set_soma_field(value, "varm", "SOMACollection")
    },

    #' @field varp a [`SOMACollection`] of [`SOMASparseNDArray`]s containing
    #' pairwise annotations on the variable axis and indexed with `[varid_1,
    #' varid_2]`.
    varp = function(value) {
      private$get_or_set_soma_field(value, "varp", "SOMACollection")
    }
  )
)
