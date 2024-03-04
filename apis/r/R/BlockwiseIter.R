#' SOMA Blockwise Read Iterator Base Class
#'
#' @description Class that allows for blockwise read iteration of SOMA reads
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
    #' @template param-blockwise-iter
    #' @template param-coords-iter
    #' @template param-dots-ignored
    #'
    initialize = function(
      sr,
      array,
      coords,
      axis,
      ...,
      reindex_disable_on_axis = NULL
    ) {
      super$initialize(sr)
      stopifnot(
        "'array' must be a 'SOMASparseNDArray'" = inherits(array, "SOMASparseNDArray")
      )
      private$.array <- array
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
        "'coords' must be a named list of 'CoordsStrider' objects" = is_named_list(coords) &&
          all(vapply_lgl(coords, inherits, "CoordsStrider"))
      )
      axname <- self$array$dimnames()[self$axis + 1L]
      if (!axname %in% names(coords)) {
        stop(
          "'coords' must include an entry for ",
          sQuote(axname, FALSE),
          call. = FALSE
        )
      }
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
    read_complete = function() !self$coords_axis$has_next() ||
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
      private$reset()
      dimnam <- self$array$dimnames()[self$axis + 1L]
      nextelems <- self$coords_axis$next_element()
      private$set_dim_points(dimnam, nextelems)
      return(private$.read_next())
    }
  ),
  active = list(
    #' @field array The underlying SOMA array
    #'
    array = function() private$.array,
    #' @field axis The axis to iterate over in a blockwise fashion
    #'
    axis = function() private$.axis,
    #' @field coords A list of \code{\link{CoordsStrider}} objects
    #'
    coords = function() private$.coords,
    #' @field coords_axis The \code{\link{CoordsStrider}} for \code{axis}
    #'
    coords_axis = function() {
      dname <- self$array$dimnames()[self$axis + 1L]
      return(self$coords[[dname]])
    },
    #' @field reindex_disable_on_axis Additional axes that will not be re-indexed
    #'
    reindex_disable_on_axis = function() private$.reindex_disable_on_axis
  ),
  private = list(
    .array = NULL,
    .coords = list(),
    .axis = integer(1L),
    .reindex_disable_on_axis = NULL,
    # @description Reset internal state of SOMA Reader while keeping array open
    reset = function() {
      if (is.null(private$soma_reader_pointer)) {
        return(NULL)
      }
      sr_reset(private$soma_reader_pointer)
      return(invisible(NULL))
    },
    # @description Set dimension selection on given axis
    set_dim_points = function(dimname, points) {
      stopifnot(
        "Name of dimension must be character" = is.character(dimname),
        "Points must be int64 vector" = inherits(points, "integer64")
      )
      if (is.null(private$soma_reader_pointer)) {
        return(NULL)
      }
      sr_set_dim_points(private$soma_reader_pointer, dimname, points)
      return(invisible(NULL))
    }
  )
)

#' SOMA Blockwise Read Iterator for Arrow Tables
#'
#' @description Class that allows for blockwise read iteration of SOMA reads
#' as Arrow \code{\link[Arrow]{Table}s}
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
#' @description Class that allows for blockwise read iteration of SOMA reads
#' as sparse matrices
#'
#' @keywords internal
#'
#' @export
#'
BlockwiseSparseReadIter <- R6::R6Class(
  classname = "BlockwiseSparseReadIter",
  inherit = BlockwiseReadIterBase,
  public = list(
    #' @description Create
    #'
    #' @template param-blockwise-iter
    #' @template param-coords-iter
    #' @template param-dots-ignored
    #' @template param-repr-read
    #'
    initialize = function(
      sr,
      array,
      coords,
      axis,
      ...,
      repr = "T",
      reindex_disable_on_axis = NULL
    ) {
      super$initialize(
        sr,
        array,
        coords,
        axis,
        ...,
        reindex_disable_on_axis = reindex_disable_on_axis
      )
      private$.repr <- match.arg(repr)
      private$.shape <- sapply(coords, length)
    },
    #' @description ...
    #'
    #' @return ...
    #'
    concat = function() soma_array_to_sparse_matrix_concat(self, private$.zero_based)
  ),
  active = list(
    #' @field repr Representation of the sparse matrix to return
    #'
    repr = function() private$.repr
  ),
  private = list(
    .repr = character(1L),
    .shape = NULL,
    .zero_based = FALSE,
    soma_reader_transform = function(x) arrow_table_to_sparse(
      soma_array_to_arrow_table(x),
      repr = self$repr,
      shape = private$.shape,
      zero_based = private$.zero_based
    )
  )
)
