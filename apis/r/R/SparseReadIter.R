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
    #' @param uri Character value with URI path to a SOMADataFrame or SOMASparseNDArray
    #' @param config character vector containing TileDB config.
    #' @param colnames Optional vector of character value with the name of the columns to retrieve
    #' @param qc Optional external Pointer object to TileDB Query Condition, defaults to \sQuote{NULL} i.e.
    #' no query condition
    #' @param dim_points Optional named list with vector of data points to select on the given
    #' dimension(s). Each dimension can be one entry in the list.
    #' @param loglevel Character value with the desired logging level, defaults to \sQuote{auto}
    #' @param repr Optional one-character code for sparse matrix representation type
    #' which lets prior setting prevail, any other value is set as new logging level.
    #' @param shape Numerical vector with two elements.
    initialize = function(uri, 
                          config, 
                          colnames = NULL, 
                          qc = NULL, 
                          dim_points = NULL, 
                          loglevel = "auto", 
                          repr = c("C", "T", "R"),
                          shape) {
        
        # Initiate super class
          super$initialize (uri = uri, config = config, colnames = colnames, qc = qc,
                            dim_points = dim_points, loglevel = loglevel)
    
          private$repr <- repr
          private$shape <- shape 
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
        mat <- mat + self$read_next()
      }
      
      mat
      
    }),
  
  private = list(
                 
    repr=NULL,
    shape=NULL,
    
    ## refined from base class
    soma_reader_transform = function(x) {
      arrow_table_to_sparse(soma_array_to_arrow_table(x), repr = private$repr, shape = private$shape)
   }
    
  )
)
