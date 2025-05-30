#' SOMA Read Iterator Over Sparse Matrices
#'
#' @description \code{SparseReadIter} is a class that allows for iteration over
#' a reads on \link{SOMASparseNDArray}.
#'
#' @export
#'
#' @examplesIf requireNamespace("withr", quietly = TRUE)
#' dir <- withr::local_tempfile(pattern = "matrix-iter")
#' dir.create(dir, recursive = TRUE)
#' (exp <- load_dataset("soma-exp-pbmc-small", dir))
#' qry <- exp$axis_query("RNA")
#' xqry <- qry$X("data")
#'
#' iter <- xqry$sparse_matrix()
#' stopifnot(inherits(iter, "SparseReadIter"))
#'
#' while (!iter$read_complete()) {
#'   block <- iter$read_next()
#' }
#'
#' \dontshow{
#' exp$close()
#' }
#'
SparseReadIter <- R6::R6Class(
  classname = "SparseReadIter",
  inherit = ReadIter,
  public = list(

    #' @description Create (lifecycle: maturing).
    #'
    #' @param sr Soma reader pointer.
    #' @param shape Shape of the full matrix.
    #' @param zero_based Logical, if \code{TRUE} will make iterator for
    #' Matrix::\link[Matrix]{dgTMatrix-class}
    #' otherwise \link{matrixZeroBasedView}.
    #'
    initialize = function(sr, shape, zero_based = FALSE) {
      # TODO implement zero_based argument, currently doesn't do anything
      stopifnot(
        "'shape' must have two dimensions" = length(shape) == 2,
        "'shape' must not exceed '.Machine$integer.max'" =
          all(shape <= .Machine$integer.max)
      )

      # Initiate super class
      super$initialize(sr)
      private$repr <- "T"
      private$shape <- shape
      private$zero_based <- zero_based
    },

    #' @description  Concatenate remainder of iterator.
    #'
    #' @return \link{matrixZeroBasedView} of Matrix::\link[Matrix]{sparseMatrix}.
    #'
    concat = function() soma_array_to_sparse_matrix_concat(
      self,
      private$zero_based
    )
  ),
  private = list(
    repr = NULL,
    shape = NULL,
    zero_based = NULL,

    # refined from base class
    #
    soma_reader_transform = function(x) {
      return(arrow_table_to_sparse(
        soma_array_to_arrow_table(x),
        repr = private$repr,
        shape = private$shape,
        zero_based = private$zero_based
      ))
    }
  )
)
