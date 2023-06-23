#' SOMASparseNDArrayRead
#'
#' @description
#' Intermediate type to choose result format when reading a sparse array
#' @export

SOMASparseNDArrayRead <- R6::R6Class(
  classname = "SOMASparseNDArrayRead",

  public = list(

    #' @description Create (lifecycle: experimental)
    #' @param sr soma read pointer
    #' @param shape Shape of the full matrix
    initialize = function(sr, shape) {
      private$sr <- sr
      private$shape <- shape
    },

    #' @description Read as a sparse matrix (lifecycle: experimental). Returns
    #' an iterator of Matrix::\link[Matrix]{dgTMatrix-class} or \link{matrixZeroBasedView} of it.
    #' @param zero_based Logical, if \code{TRUE} returns iterator of \link{matrixZeroBasedView}
    #' if \code{FALSE} returns iterator of Matrix::\link[Matrix]{dgTMatrix-class}.
    #' @return \link{SparseReadIter}
    sparse_matrix = function(zero_based=FALSE) {
      #TODO implement zero_based argument, currently doesn't do anything


      if (any(private$shape >= .Machine$integer.max)) {
        warning(
          "Array's 0-based domain exceeds '.Machine$integer.max'.\n",
          "  - Result will only include coordinates within [0, 2^31 - 1).\n",
          "  - The full range of coordinates can be obtained with $tables().",
          call. = FALSE,
          immediate. = TRUE
        )
        private$shape <- pmin(private$shape, .Machine$integer.max - 1L)
      }

      SparseReadIter$new(private$sr, private$shape, zero_based = zero_based)
    },

    #' @description Read as a arrow::\link[arrow]{Table} (lifecycle: experimental).
    #' Returns an iterator of arrow::\link[arrow]{Table}.
    #' @return \link{TableReadIter}
    tables = function() {
      TableReadIter$new(private$sr)
    }
  ),

  private = list(
    sr=NULL,
    shape=NULL
  )

)
