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
    concat = function() soma_array_to_arrow_table_concat(self)
  ),

  private = list(
    ## refined from base class
    soma_reader_transform = function(x) {
      at <- soma_array_to_arrow_table(x)
      at
    }
  )
)
