#' SOMA Read Iterator Over Arrow Tables
#'
#' @description \code{TableReadIter} is a class that allows for iteration over
#' a reads on \link{SOMASparseNDArray} and \link{SOMADataFrame}.
#' Iteration chunks are retrieved as an Arrow \link[arrow]{Table}
#'
#' @export
#'
#' @examplesIf requireNamespace("withr", quietly = TRUE)
#' dir <- withr::local_tempfile(pattern = "blockwise-table")
#' dir.create(dir, recursive = TRUE)
#' (exp <- load_dataset("soma-exp-pbmc-small", dir))
#' qry <- exp$axis_query("RNA")
#' xqry <- qry$X("data")
#'
#' iter <- xqry$tables()
#' stopifnot(inherits(iter, "TableReadIter"))
#'
#' while (!iter$read_complete()) {
#'   block <- iter$read_next()
#' }
#'
#' \dontshow{
#' exp$close()
#' }
#'
TableReadIter <- R6::R6Class(
  classname = "TableReadIter",
  inherit = ReadIter,
  public = list(
    #' @description  Concatenate remainder of iterator
    #'
    #' @return An Arrow \code{\link[arrow:Table]{Table}}
    #'
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
