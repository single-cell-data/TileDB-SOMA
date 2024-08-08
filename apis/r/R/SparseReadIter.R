#' SparseReadIter
#'
#' @description
#' \code{SparseReadIter} is a class that allows for iteration over
#'  a reads on \link{SOMASparseNDArray}.
#' Iteration chunks are retrieved as 0-based Views \link{matrixZeroBasedView} of Matrix::\link[Matrix]{sparseMatrix}.
#' @export

SparseReadIter <- R6::R6Class(
  classname = "SparseReadIter",
  inherit = ReadIter,

  public = list(

    #' @description Create (lifecycle: maturing)
    #' @param sr Soma reader pointer
    #' @param shape Shape of the full matrix
    #' @param zero_based Logical, if TRUE will make iterator for Matrix::\link[Matrix]{dgTMatrix-class}
    #' otherwise \link{matrixZeroBasedView}.
    initialize = function(sr, shape, zero_based=FALSE) {
      #TODO implement zero_based argument, currently doesn't do anything
      stopifnot("'shape' must have two dimensions" = length(shape) == 2,
                "'shape' must not exceed '.Machine$integer.max'" =
                  all(shape <= .Machine$integer.max))

      # Initiate super class
      super$initialize(sr)
      private$repr <- "T"
      private$shape <- shape
      private$zero_based <- zero_based
    },

    #' @description  Concatenate remainder of iterator.
    #' @return \link{matrixZeroBasedView} of Matrix::\link[Matrix]{sparseMatrix}
    concat = function() soma_array_to_sparse_matrix_concat(self, private$zero_based)
    ),

  private = list(

    repr=NULL,
    shape=NULL,
    zero_based=NULL,

    ## refined from base class
    soma_reader_transform = function(x) {
      arrow_table_to_sparse(soma_array_to_arrow_table(x),
                            repr = private$repr,
                            shape = private$shape,
                            zero_based = private$zero_based)
   }

  )
)
