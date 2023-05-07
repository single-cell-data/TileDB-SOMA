#' TableReadFull
#'
#' @description TODO
#' @export

TableReadFull <- R6::R6Class(
  classname = "TableReadFull",
  inherit = ReadFull,
  
  private = list(
    ## refined from base class
    soma_reader_transform = function(x) {
      arrow::as_arrow_table(arrow::RecordBatch$import_from_c(x[[1]], x[[2]]))
    }
    
  )
)
