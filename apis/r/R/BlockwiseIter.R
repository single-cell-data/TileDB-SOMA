#' SOMA Blockwise Read Iterator Base Class
#'
#' #' @description Class that allows for blockwise read iteration of SOMA reads
#'
#' @keywords internal
#'
#' @export
#'
BlockwiseReadIterBase <- R6::R6Class(
  classname = "BlockwiseReadIterBase",
  inherit = ReadIter,
  public = list(
    #' @description Create
    #'
    #' @param sr SOMA read pointer
    #' @param array ...
    #' @param coords ...
    #' @param axis ...
    #' @param ... Ignored
    #' @param reindex_disable_on_axis ...
    #' @param eager ...
    #'
    initialize = function(
      sr,
      array,
      coords,
      axis,
      ...,
      reindex_disable_on_axis = NULL,
      eager = TRUE
    ) {
      super$initialize(sr)
      stopifnot(
        "'array' must be a 'SOMASparseNDArray'" = inherits(array, "SOMASparseNDArray"),
        "'eager' must be TRUE or FALSE" = isTRUE(eager) || isFALSE(eager)
      )
      private$.array <- array
      private$.eager <- eager
      # Check axis
      ndim <- self$array$ndim() - 1L
      if (!rlang::is_integerish(axis, 1L, TRUE) && axis >= 0L && axis <= ndim) {
        stop(
          "'axis' must be a single integer value between 0 and ",
          ndim,
          ", inclusive",
          call. = FALSE
        )
      }
      private$.axis <- axis
      # Check coords
      stopifnot(
        is_named_list(coords) && all(vapply_lgl(coords, inherits, "CoordsStrider")),
        self$array$dimnames()[self$axis + 1L] %in% names(coords)
      )
      private$.coords <- coords
      # Check reindex_disable_on_axis
      if (!is.null(reindex_disable_on_axis)) {
        stopifnot(
          "'reindex_disable_on_axis' must be a vector of integers" = (
            rlang::is_integerish(reindex_disable_on_axis) ||
              inherits(reindex_disable_on_axis, "integer64")
          ),
          "'reindex_disable_on_axis' must be finite" = is.finite(reindex_disable_on_axis),
          "'reindex_disable_on_axis' must be within the range of dimensions of the array" = all(
            reindex_disable_on_axis >= 0 && reindex_disable_on_axis <= ndim
          )
        )
      }
      private$.reindex_disable_on_axis <- reindex_disable_on_axis
    },
    #' @description Check if the iterated read is complete or not
    #'
    #' @return \code{TRUE} if read is complete, otherwise \code{FALSE}
    #'
    read_complete = function() !self$coords_axis$hasNext() ||
      is.null(private$soma_reader_pointer),
    #' @description Read the next chunk of the iterated read. If read
    #' is complete, throws an \code{iterationCompleteWarning} warning and
    #' returns \code{NULL}
    #'
    #' @return \code{NULL} or the next blockwise chunk of the iterated read
    #'
    read_next = function() {
      if (is.null(private$soma_reader_pointer)) {
       return(NULL)
      }
      if (self$read_complete()) {
        return(private$.readComplete())
      }
      super$reset()
      dimnam <- self$array$dimnames()[self$axis + 1L]
      nextelems <- self$coords_axis$nextElem()
      super$set_dim_points(dimnam, nextelems)
      val <- private$.read_next()
      val
    }
  ),
  active = list(
    #' @field array The underlying SOMA array
    array = function() private$.array,
    #' @field context Not yet implemented
    context = function() .NotYetImplemented(),
    #' @field axis The axis to iterate over in a blockwise fashion
    axis = function() private$.axis,
    #' @field coords A list of \code{\link{CoordsStrider}} objects
    coords = function() private$.coords,
    #' @field coords_axis The \code{\link{CoordsStrider}} for \code{axis}
    coords_axis = function() {
      dname <- self$array$dimnames()[self$axis + 1L]
      return(self$coords[[dname]])
    },
    #' @field reindex_disable_on_axis ...
    reindex_disable_on_axis = function() private$.reindex_disable_on_axis,
    #' @field eager ...
    eager = function() private$.eager
  ),
  private = list(
    .array = NULL,
    .coords = list(),
    .axis = integer(1L),
    .reindex_disable_on_axis = NULL,
    .eager = logical(1L)
  )
)

#' SOMA Blockwise Read Iterator for Arrow Tables
#'
#' @description ...
#'
#' @keywords internal
#'
#' @export
#'
BlockwiseTableReadIter <- R6::R6Class(
  classname = "BlockwiseTableReadIter",
  inherit = BlockwiseReadIterBase,
  public = list(
    #' @description ...
    #'
    #' @return ...
    #'
    concat = function() soma_array_to_arrow_table_concat(self)
  ),
  private = list(
    soma_reader_transform = function(x) soma_array_to_arrow_table(x)
  )
)

#' SOMA Blockwise Read Iterator for Sparse Matrices
#'
#' @description ...
#'
#' @keywords internal
#'
#' @export
#'
BlockwiseSparseReadIter <- R6::R6Class(
  classname = "BlockwiseSparseReadIter",
  inherit = BlockwiseReadIterBase,
  public = list(
    #' @description ...
    #'
    #' @param sr ...
    #' @param array ...
    #' @param coords ...
    #' @param axis ...
    #' @param ... Ignored
    #' @param repr ...
    #' @param compress ...
    #' @param reindex_disable_on_axis ...
    #' @param eager ...
    #'
    initialize = function(
      sr,
      array,
      coords,
      axis,
      ...,
      repr = "T",
      compress = TRUE,
      reindex_disable_on_axis = NULL,
      eager = TRUE
    ) {
      super$initialize(
        sr,
        array,
        coords,
        axis,
        ...,
        reindex_disable_on_axis = reindex_disable_on_axis,
        eager = eager
      )
      private$.repr <- match.arg(repr)
      stopifnot(isTRUE(compress) || isFALSE(compress))
      private$.compress <- compress
    },
    #' @description ...
    #'
    #' @return ...
    #'
    concat = function() {
      # TODO: implement concat() of blockwise sparse matrix iterator
      .NotYetImplemented()
    }
  ),
  active = list(
    #' @field repr ...
    #'
    repr = function() private$.repr,
    #' @field compress ...
    #'
    compress = function() private$.compress
  ),
  private = list(
    .repr = character(1L),
    .compress = logical(1L),
    soma_reader_transform = function(x) {
      # TODO: implement soma_reader_transform() of blockwise sparse matrix iterator
      .NotYetImplemented()
    }
  )
)
