#' SparseReadIter
#'
#' @description
#' \code{SparseReadIter} is a class that allows for iteration over 
#'  a reads on \link{SOMASparseNDArray}.
#' Iteration chunks are retrieved as 0-based Views \link{matrixZeroBasedView} of Matrix::\link[Matrix]{SparseMatrix}.
#' @export

SparseReadIter <- R6::R6Class(
  classname = "SparseReadIter",
  inherit = ReadIter,

  public = list(
                
    #' @description Create (lifecycle: experimental)
    #' @param shape Shape of the full matrix
    #' @param zero_based Logical, if TRUE will make iterator for Matrix::\link[Matrix]{dgTMatrix}
    #' otherwise \link{matrixZeroBasedView}.
    initialize = function(sr, shape, zero_based=FALSE) {
      #TODO implement zero_based argument, currently doesn't do anything
      stopifnot("Array must have two dimensions" = length(shape) == 2,
                "Array dimensions must not exceed '.Machine$integer.max'" = any(shape < .Machine$integer.max))
      
      # Initiate super class
        super$initialize(sr)
        private$repr <- "T"
        private$shape <- shape 
        private$zero_based <- zero_based
    },
    
   
    #' @description  Concatenate remainder of iterator.
    #' @return \link{matrixZeroBasedView} of Matrix::\link[Matrix]{SparseMatrix}
    concat = function(){
      
      if(self$read_complete()) {
        warning("Iteration complete, returning NULL")
        return(NULL)
      }
      
      mat <- self$read_next()
      
      while (!self$read_complete()) {
        if(private$zero_based) {
            mat <- mat$sum(self$read_next())
        } else {
            mat <- mat + self$read_next()
        }
      }
      
      mat
      
    }),
  
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
