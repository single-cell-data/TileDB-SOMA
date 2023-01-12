#' SOMA Measurement
#'
#' @description A `SOMAMeasurement` is a sub-element of a [`SOMAExperiment`], and is otherwise a specialized [`SOMACollection`] with pre-defined fields:
#'
#' - `var` ([`SOMADataFrame`]): Primary annotations on the variable axis, for
#'   variables in this measurement (i.e., annotates columns of `X`). The
#'   contents of the `soma_joinid` column define the variable index domain, AKA
#'   var_id. All variables for this measurement must be defined in this
#'   dataframe.
#' - `X` ([`SOMACollection`] of [`SOMASparseNDArray`]): A collection of sparse
#'   matrices, each containing measured feature values. Each matrix is indexed
#'   by `[obsid, varid]`.
#' - `obsm`: ([`SOMACollection`] of [`DenseNDArray`]): A collection of dense
#'   matrices containing annotations of each `obs` row. Has the same shape as
#'   `obs`, and is indexed with `obsid`.
#' - `obsp`: ([`SOMACollection`] of [`SparseNDArray`]): A collection of sparse
#'   matrices containing pairwise annotations of each `obs` row. Indexed with
#'   `[obsid_1, obsid_2]`.
#' - `varm`: ([`SOMACollection`] of [`DenseNDArray`]): A collection of dense
#'   matrices containing annotations of each `var` row. Has the same shape as
#'   `var`, and is indexed with `varid`.
#' - `varp`: ([`SOMACollection`] of [`SparseNDArray`]): A collection of sparse
#'   matrices containing pairwise annotations of each `var` row. Indexed with
#'   `[varid_1, varid_2]`
#'
#' @export
SOMAMeasurement <- R6::R6Class(
  classname = "SOMAMeasurement",
  inherit = SOMACollectionBase,

  active = list(
    #' @field Retrieve of set `var` [`SOMADataFrame`].
    var = function(value) {
      private$get_or_set_soma_field(value, "var", "SOMADataFrame")
    },

    #' @field Retrieve or set `X` [`SOMACollection`].
    X = function(value) {
      private$get_or_set_soma_field(value, "X", "SOMACollection")
    },

    #' @field Retrieve or set `obsm` [`SOMACollection`].
    obsm = function(value) {
      private$get_or_set_soma_field(value, "obsm", "SOMACollection")
    },

    #' @field Retrieve or set `obsp` [`SOMACollection`].
    obsp = function(value) {
      private$get_or_set_soma_field(value, "obsp", "SOMACollection")
    },

    #' @field Retrieve or set `varm` [`SOMACollection`].
    varm = function(value) {
      private$get_or_set_soma_field(value, "varm", "SOMACollection")
    },

    #' @field Retrieve or set `varp` [`SOMACollection`].
    varp = function(value) {
      private$get_or_set_soma_field(value, "varp", "SOMACollection")
    }
  )
)
