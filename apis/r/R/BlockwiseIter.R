BlockwiseReadIterBase <- R6::R6Class(
  classname = "BlockwiseReadIterBase",
  inherit = ReadIter,
  public = list(
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
    read_complete = function() !self$coords_axis$hasNext() || super$read_complete(),
    read_next = function() {
      message("blockwise read next")
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
      val <- super$read_next()
      #print(str(val))  # TODO
      val
    }
  ),
  active = list(
    array = function() private$.array,
    context = function() .NotYetImplemented(),
    axis = function() private$.axis,
    coords = function() private$.coords,
    coords_axis = function() {
      dname <- self$array$dimnames()[self$axis + 1L]
      return(self$coords[[dname]])
    },
    reindex_disable_on_axis = function() private$.reindex_disable_on_axis,
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

BlockwiseTableReadIter <- R6::R6Class(
  classname = "BlockwiseTableReadIter",
  inherit = BlockwiseReadIterBase,
  public = list(
    concat = function() {
      message("blockwise table concat")
      # TODO: implement concat() of blockwise table iterator
      .NotYetImplemented()
    }
  ),
  private = list(
    soma_reader_transform = function(x) {
      message("blockwise table transform")
      # TODO: implement soma_reader_transform() of blockwise table iterator
      .NotYetImplemented()
    }
  )
)

BlockwiseSparseReadIter <- R6::R6Class(
  classname = "BlockwiseSparseReadIter",
  inherit = BlockwiseReadIterBase,
  public = list(
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
    concat = function() {
      message("blockwise sparse concat")
      # TODO: implement concat() of blockwise sparse matrix iterator
      .NotYetImplemented()
    }
  ),
  active = list(
    repr = function() private$.repr,
    compress = function() private$.compress
  ),
  private = list(
    .repr = character(1L),
    .compress = logical(1L),
    soma_reader_transform = function(x) {
      message("blockwise sparse transform")
      # TODO: implement soma_reader_transform() of blockwise sparse matrix iterator
      .NotYetImplemented()
    }
  )
)
