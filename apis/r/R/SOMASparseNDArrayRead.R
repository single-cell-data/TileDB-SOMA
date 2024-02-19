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
    #' @param sr SOMA read pointer
    #' @param array \code{\link{SOMASparseNDArray}}
    #' @param coords ...
    # @param shape Shape of the full matrix
    #'
    initialize = function(sr, array, coords) {
      stopifnot(
        "'array' must be a SOMASparseNDArray" = inherits(array, "SOMASparseNDArray")
      )
      if (is.null(coords)) {
        private$.coords <- vector(mode = "list", length = array$ndim())
        shape <- array$shape()
        for (i in seq_along(private$.coords)) {
          private$.coords[[i]] <- CoordsStrider$new(
            start = 0L,
            end = shape[i],
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
          coords_vec <- sapply(
            X = coords,
            FUN = function(x) x$coords,
            simplify = FALSE,
            USE.NAMES = TRUE
          )
          Filter(Negate(is.null), x = coords_vec)
          if (length(coords_vec)) {
            private$.coords_vec <- coords_vec
          }
        } else {
          private$.coords_vec <- coords
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
    sr = function() return(private$.sr),
    #' @field array The underlying \code{\link{SOMASparseNDArray}}
    array = function() return(private$.array),
    #' @field coords The iterated coordinates for the read
    coords = function() return(private$.coords),
    #' @field coords_vec If \code{coords} is passed, then the coordinate vector
    #'
    coords_vec = function() return(private$.coords_vec),
    #' @field shape The shape of the underlying array
    shape = function() return(self$array$shape())
  ),
  private = list(
    .sr = NULL,
    .array = NULL,
    .coords = NULL,
    .coords_vec = NULL
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

    # @description Create (lifecycle: experimental)
    # @param sr soma read pointer
    # @param shape Shape of the full matrix
    # initialize = function(sr, shape) {
    #   private$sr <- sr
    #   private$shape <- shape
    # },

    #' @description Read as a sparse matrix (lifecycle: experimental). Returns
    #' an iterator of Matrix::\link[Matrix]{dgTMatrix-class} or \link{matrixZeroBasedView} of it.
    #' @param zero_based Logical, if \code{TRUE} returns iterator of \link{matrixZeroBasedView}
    #' if \code{FALSE} returns iterator of Matrix::\link[Matrix]{dgTMatrix-class}.
    #' @return \link{SparseReadIter}
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

    #' @description Read as a arrow::\link[arrow]{Table} (lifecycle: experimental).
    #' Returns an iterator of arrow::\link[arrow]{Table}.
    #' @return \link{TableReadIter}
    tables = function() {
      TableReadIter$new(self$sr)
    },
    #' @description ...
    #'
    #' @param axis ...
    #' @param ... Ignored
    #' @param size ...
    #' @param reindex_disable_on_axis ...
    #' @param eager ...
    #'
    #' @return A \code{\link{SOMASparseNDArrayBlockwiseRead}} iterated reader
    #'
    blockwise = function(
      axis,
      ...,
      size = NULL,
      reindex_disable_on_axis = NULL,
      eager = TRUE
    ) {
      return(SOMASparseNDArrayBlockwiseRead$new(
        self$sr,
        self$array,
        self$coords,
        axis,
        size = size,
        reindex_disable_on_axis = reindex_disable_on_axis,
        eager = eager
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
    #' @param sr SOMA read pointer
    #' @param array \code{\link{SOMASparseNDArray}}
    #' @param coords ...
    #' @param axis ...
    #' @param ... Ignored
    #' @param size ...
    #' @param reindex_disable_on_axis ...
    #' @param eager ...
    #'
    initialize = function(
      sr,
      array,
      coords,
      axis,
      ...,
      size,
      reindex_disable_on_axis,
      eager = TRUE
    ) {
      super$initialize(sr, array, coords)
      stopifnot(
        is.null(size) ||
          rlang::is_integerish(size, finite = TRUE) ||
          (inherits(size, 'integer64') && all(is.finite(size))),
        is.null(reindex_disable_on_axis) ||
          rlang::is_integerish(reindex_disable_on_axis, finite = TRUE) ||
          (inherits(reindex_disable_on_axis, 'integer64') && all(is.finite(reindex_disable_on_axis))),
        isTRUE(eager) || isFALSE(eager)
      )
      private$.axis <- axis
      private$.size <- size
      private$.reindex_disable_on_axis <- reindex_disable_on_axis
      private$.eager <- eager
    },
    #' @description ...
    #'
    #' @return ...
    #'
    tables = function() {
      .NotYetImplemented()
    },
    #' @description ...
    #'
    #' @param compress ...
    #'
    #' @return ...
    #'
    sparse_matrix = function(compress = TRUE) {
      stopifnot(
        "'compress' must be TRUE or FALSE" = isTRUE(compress) || isFALSE(compress)
      )
      .NotYetImplemented()
    }
  ),
  active = list(
    #' @field axis ...
    axis = function() return(private$.axis),
    #' @field size ...
    size = function() return(private$.size),
    #' @field reindex_disable_on_axis ...
    reindex_disable_on_axis = function() return(private$.reindex_disable_on_axis),
    #' @field eager ...
    eager = function() return(private$eager)
  ),
  private = list(
    .axis = NULL,
    .size = NULL,
    .reindex_disable_on_axis = NULL,
    .eager = NULL
  )
)
