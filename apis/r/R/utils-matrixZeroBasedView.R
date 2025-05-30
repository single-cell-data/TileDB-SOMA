#' Zero-based Wrapper for Sparse Matrices
#'
#' \code{matrixZeroBasedView} is a wrapper shim for a matrix or
#' \code{\link[Matrix:sparseMatrix]{Matrix::sparseMatrix}} that allows
#' elemental matrix access using zero-based indices.
#'
#' @export
#'
#' @examples
#' (mat <- Matrix::rsparsematrix(3L, 3L, 0.3))
#' (mat0 <- matrixZeroBasedView$new(mat))
#'
#' mat0$take(0, 0)
#' mat0$take(0, 0:2)$get_one_based_matrix()
#'
matrixZeroBasedView <- R6::R6Class(
  classname = "matrixZeroBasedView",
  public = list(
    #' @description Initialize (lifecycle: maturing).
    #'
    #' @param x A matrix.
    #'
    initialize = function(x) {
      if (!inherits(x, "matrix") && !inherits(x, "sparseMatrix") && !inherits(one_based_matrix, "Matrix")) {
        stop("Matrix object must inherit class matrix or Matrix::sparseMatrix or Matrix::Matrix")
      }
      if (length(dim(x)) != 2) {
        stop("Only two-dimensional matrices are supported")
      }
      private$one_based_matrix <- x
    },

    #' @description Zero-based matrix element access.
    #'
    #' @param i Row index (zero-based).
    #' @param j Column index (zero-based).
    #'
    #' @return The specified matrix slice as another \code{matrixZeroBasedView}.
    #'
    take = function(i = NULL, j = NULL) {
      x <- NULL
      if (is.null(i) && is.null(j)) {
        x <- private$one_based_matrix[, , drop = FALSE]
      } else if (is.null(i)) {
        x <- private$one_based_matrix[, j + 1, drop = FALSE]
      } else if (is.null(j)) {
        x <- private$one_based_matrix[i + 1, , drop = FALSE]
      } else {
        x <- private$one_based_matrix[i + 1, j + 1, drop = FALSE]
      }

      matrixZeroBasedView$new(x)
    },

    #' @description dim.
    #'
    #' @return The dimensions of the matrix.
    #'
    dim = function() {
      dim(private$one_based_matrix)
    },

    #' @description nrow.
    #'
    #' @return Matrix row count.
    #'
    #' @export
    #'
    nrow = function() {
      nrow(private$one_based_matrix)
    },

    #' @description ncol.
    #'
    #' @return Matrix column count.
    #'
    ncol = function() {
      ncol(private$one_based_matrix)
    },

    #' @description Get the one-based \R matrix with its original class.
    #'
    #' @return One-based matrix.
    #'
    get_one_based_matrix = function() {
      private$one_based_matrix
    },

    #' @description Perform arithmetic sum between this
    #' \code{matrixZeroBasedView} and another \code{matrixZeroBasedView}.
    #'
    #' @param x the \code{matrixZeroBasedView} to sum.
    #'
    #' @return The result of the sum as a \code{matrixZeroBasedView}.
    #'
    sum = function(x) {
      if (!inherits(x, "matrixZeroBasedView")) {
        stop("Only arithmetic sum with another 'matrixZeroBasedView` is supported")
      }
      matrixZeroBasedView$new(private$one_based_matrix + x$get_one_based_matrix())
    },

    #' @description print.
    #'
    #' @return Invisibly returns \code{self}.
    #'
    print = function() {
      dims <- self$dim()
      cat("Non-mutable 0-based 'view' class for matrices.\n")
      cat("To get 1-based matrix use `x$get_one_based_matrix()\n")
      cat(paste0("Dimensions: ", dims[1], "x", dims[2], "\n"))
      return(invisible(self))
    }
  ),
  private = list(
    one_based_matrix = NULL
  )
)
