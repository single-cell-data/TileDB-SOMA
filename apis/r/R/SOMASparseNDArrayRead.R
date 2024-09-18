#' SOMA Sparse ND-Array Reader Base
#'
#' @description Base class for SOMA sparse ND-array reads
#'
#' @keywords internal
#'
#' @export
#'
SOMASparseNDArrayReadBase <- R6::R6Class(
  classname = "SOMASparseNDArrayReadBase",
  cloneable = FALSE,
  public = list(
    #' @description Create
    #'
    #' @template param-blockwise-iter
    #' @template param-coords-read
    #'
    initialize = function(sr, array, coords = NULL) {
      stopifnot(
        "'array' must be a SOMASparseNDArray" = inherits(array, "SOMASparseNDArray")
      )
      if (is.null(coords)) {
        private$.coords <- vector(mode = "list", length = array$ndim())
        shape <- array$shape()
        for (i in seq_along(private$.coords)) {
          private$.coords[[i]] <- CoordsStrider$new(
            start = 0L,
            end = shape[i] - 1L,
            stride = .Machine$integer.max
          )
        }
        names(private$.coords) <- array$dimnames()
      } else {
        stopifnot(
          "'coords' must be a list of integer64 values" = is.list(coords) &&
            all(vapply_lgl(coords, inherits, what = c('integer64', 'numeric', 'CoordsStrider'))),
          "'coords' must be named with the dimnames of 'array'" = is_named(coords, FALSE) &&
            all(names(coords) %in% array$dimnames())
        )
        if (all(vapply_lgl(coords, inherits, what = 'CoordsStrider'))) {
          private$.coords <- coords
        } else {
          private$.coords <- vector(mode = "list", length = length(coords))
          names(private$.coords) <- names(coords)
          for (i in names(coords)) {
            private$.coords[[i]] <- CoordsStrider$new(coords[[i]], stride = .Machine$integer.max)
          }
        }
      }
      private$.sr <- sr
      private$.array <- array
    }
  ),
  active = list(
    #' @field sr The SOMA read pointer
    #'
    sr = function() return(private$.sr),
    #' @field array The underlying \code{\link{SOMASparseNDArray}}
    #'
    array = function() return(private$.array),
    #' @field coords The iterated coordinates for the read
    #'
    coords = function(value) {
      if (missing(value)) {
        return(private$.coords)
      }
      if (!is.list(value) && is_named(value, allow_empty = FALSE)) {
        stop("'coords' must be a named list", call. = FALSE)
      }
      if (!all(names(x = value) %in% names(private$.coords))) {
        stop(
          "'coords' must be named with ",
          paste(sQuote(names(private$.coords)), collapse = ', '),
          call. = FALSE
        )
      }
      if (!all(vapply_lgl(value, inherits, what = 'CoordsStrider'))) {
        stop("'coords' must be a list of CoordsStriders", call. = FALSE)
      }
      for (dim in names(value)) {
        strider <- value[[dim]]
        current <- private$.coords[[dim]]
        checks <- c(
          start = strider$start == current$start,
          end = strider$end == current$end,
          coords = identical(strider$coords, current$coords)
        )
        if (!all(checks)) {
          stop(
            "New striders must cover the same coordinates as existing striders (offending: ",
            sQuote(dim),
            ")",
            call. = FALSE
          )
        }
        private$.coords[[dim]] <- strider
      }
      return(invisible(NULL))
    },
    #' @field shape The shape of the underlying array
    #'
    shape = function() return(self$array$shape())
  ),
  private = list(
    .sr = NULL,
    .array = NULL,
    .coords = NULL
  )
)

#' SOMASparseNDArrayRead
#'
#' @description
#' Intermediate type to choose result format when reading a sparse array
#' @keywords internal
#' @export

SOMASparseNDArrayRead <- R6::R6Class(
  classname = "SOMASparseNDArrayRead",
  inherit = SOMASparseNDArrayReadBase,
  cloneable = FALSE,
  public = list(

    #' @description Read as a sparse matrix (lifecycle: maturing). Returns
    #' an iterator of Matrix::\link[Matrix]{dgTMatrix-class} or
    #' \link{matrixZeroBasedView} of it.
    #'
    #' @param zero_based Logical, if \code{TRUE} returns iterator of
    #' \link{matrixZeroBasedView} if \code{FALSE} returns iterator of
    #' Matrix::\link[Matrix]{dgTMatrix-class}.
    #'
    #' @return \link{SparseReadIter}
    #'
    sparse_matrix = function(zero_based=FALSE) {
      #TODO implement zero_based argument, currently doesn't do anything

      shape <- self$shape
      # if (any(private$shape > .Machine$integer.max)) {
      if (any(shape > .Machine$integer.max)) {
        warning(
          "Array's shape exceeds '.Machine$integer.max'.\n",
          "  - Result will only include coordinates within [0, 2^31 - 1).\n",
          "  - The full range of coordinates can be obtained with $tables().",
          call. = FALSE,
          immediate. = TRUE
        )
        # private$shape <- pmin(private$shape, .Machine$integer.max)
        shape <- pmin(shape, .Machine$integer.max)
      }
      SparseReadIter$new(self$sr, shape, zero_based = zero_based)
    },

    #' @description Read as a arrow::\link[arrow]{Table} (lifecycle: maturing).
    #' Returns an iterator of arrow::\link[arrow]{Table}.
    #'
    #' @return \link{TableReadIter}
    #'
    tables = function() {
      TableReadIter$new(self$sr)
    },
    #' @description Read in a blockwise fashion
    #'
    #' @template param-blockwise-iter
    #' @template param-dots-ignored
    #'
    #' @return A \code{\link{SOMASparseNDArrayBlockwiseRead}} iterated reader
    #'
    blockwise = function(
      axis,
      ...,
      size = NULL,
      reindex_disable_on_axis = NA
    ) {
      return(SOMASparseNDArrayBlockwiseRead$new(
        self$sr,
        self$array,
        self$coords,
        axis,
        size = size,
        reindex_disable_on_axis = reindex_disable_on_axis
      ))
    }
  )
)

#' Blockwise Sparse ND-Array Reader
#'
#' @description Blockwise reader for \code{\link{SOMASparseNDArray}}
#'
#' @keywords internal
#'
#' @export
#'
SOMASparseNDArrayBlockwiseRead <- R6::R6Class(
  classname = "SOMASparseNDArrayBlockwiseRead",
  inherit = SOMASparseNDArrayReadBase,
  cloneable = FALSE,
  public = list(
    #' @description Create
    #'
    #' @template param-blockwise-iter
    #' @template param-coords-read
    #' @template param-dots-ignored
    #'
    initialize = function(
      sr,
      array,
      coords,
      axis,
      ...,
      size,
      reindex_disable_on_axis = NA
    ) {
      super$initialize(sr, array, coords)
      stopifnot(
        "'axis' must be a single integer value" = rlang::is_integerish(axis, n = 1L, finite = TRUE),
        "'size' must be a single integer value" = is.null(size) ||
          rlang::is_integerish(size, 1L, finite = TRUE) ||
          (inherits(size, 'integer64') && length(size) == 1L && is.finite(size)),
        "'reindex_disable_on_axis' must be a vector of integers" = is.null(reindex_disable_on_axis) ||
          is_scalar_logical(reindex_disable_on_axis) ||
          rlang::is_integerish(reindex_disable_on_axis, finite = TRUE) ||
          (inherits(reindex_disable_on_axis, 'integer64') && all(is.finite(reindex_disable_on_axis)))
      )
      if (axis < 0L || axis >= self$array$ndim()) {
        stop(
          "'axis' must be between 0 and ",
          self$array$ndim() - 1L,
          call. = FALSE
        )
      }
      private$.axis <- axis
      if (!is.null(size)) {
        for (i in seq_along(self$coords)) {
          self$coords[[i]]$stride <- size
        }
      }
      private$.reindex_disable_on_axis <- reindex_disable_on_axis
    },
    #' @description Read as an \code{\link[Arrow:Table]{Arrow::Table}}
    #'
    #' @return A blockwise iterator yielding chunks as
    #' \code{\link[Arrow:Table]{Arrow::Table}s}
    #'
    tables = function() {
      return(BlockwiseTableReadIter$new(
        sr = self$sr,
        array = self$array,
        coords = self$coords,
        axis = self$axis,
        reindex_disable_on_axis = private$.reindex_disable_on_axis
      ))
    },
    #' @description Read as a sparse matrix
    #'
    #' @template param-repr-read
    #'
    #' @return A blockwise iterator yielding chunks as sparse matrices
    #'
    sparse_matrix = function(repr = "T") {
      return(BlockwiseSparseReadIter$new(
        sr = self$sr,
        array = self$array,
        coords = self$coords,
        axis = self$axis,
        repr = repr,
        reindex_disable_on_axis = private$.reindex_disable_on_axis
      ))
    }
  ),
  active = list(
    #' @field axis The axis to iterate over in a blockwise fashion
    axis = function() private$.axis
  ),
  private = list(
    .axis = NULL,
    .reindex_disable_on_axis = NULL
  )
)
