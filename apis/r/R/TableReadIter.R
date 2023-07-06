#' SOMA Read Iterator over Arrow Table
#'
#' @description
#' `TableReadIter` is a class that allows for iteration over 
#'  a reads on \link{SOMASparseNDArray} and \link{SOMADataFrame}. 
#' Iteration chunks are retrieved as arrow::\link[arrow]{Table}
#' @export

TableReadIter <- R6::R6Class(
  classname = "TableReadIter",
  inherit = ReadIter,

  public = list(
    
    #' @description  Concatenate remainder of iterator.
    #' @return arrow::\link[arrow]{Table}
    concat = function(){
  
      if(self$read_complete()) {
        warning("Iteration complete, returning NULL")
        return(NULL)
      }
    
      tbl <- self$read_next()
      
      while (!self$read_complete()) {
        tbl <- arrow::concat_tables(tbl, self$read_next())
      }
      
      tbl
      
    }),
  
  private = list(
                 
    ## refined from base class
    soma_reader_transform = function(x) {
      at <- soma_array_to_arrow_table(x)
      for (n in names(private$enums)) {     # DE::TEMP::FIXME start
          if (!is.null(private$enums[[n]])) {
              at[[n]] <- arrow::DictionaryArray$create(at[[n]]$as_vector(), private$enums[[n]])
          }
      }  # DE::TEMP::FIXME start
      at
    }

  )
)
