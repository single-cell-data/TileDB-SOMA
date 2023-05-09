#' SparseReadIter
#'
#' @description
#' `SparseReadIter` is a class that allows for iteration over 
#'  the results of a read operation from SOMA objects#' 
#' @export

SparseReadIter <- R6::R6Class(
  classname = "SparseReadIter",
  inherit = ReadIter,

  public = list(
    
    ## refined from base class
    concat = function(){
  
      rl <- list()
      
      while (!self$read_complete()) {
        rl <- c(rl, self$read_next())
      }
      
      do.call(arrow::concat_tables, rl)
      
    }),
  
  private = list(
                 
    ## refined from base class
    soma_reader_transform = function(x) {
      soma_array_to_arrow(x)
    }
    
  )
)
