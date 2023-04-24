#' TableReadIter
#'
#' @description
#' `TableReadIter` is a class that allows for iteration over 
#  the results of a read operation from SOMA objects#' @importFrom stats setNames
#' @export

TableReadIter <- R6::R6Class(
  classname = "TableReadIter",
  inherit = ReadIter,

  public = list(
    
    ## refined from base class
    read_concat = function(){
  
      rl <- list()
      
      while (!self$read_complete()) {
        rl <- c(rl, self$read_next())
      }
      
      do.call(arrow::concat_tables, rl)
      
    }),
  
  private = list(
                 
    ## refined from base class
    soma_reader_transform = function(x) {
      arrow::as_arrow_table(arrow::RecordBatch$import_from_c(x[[1]], x[[2]]))
    }
    
  )
)
