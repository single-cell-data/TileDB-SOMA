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
      reindex_disable_on_axis = NA
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
      if (is_scalar_logical(reindex_disable_on_axis)) {
        reindex_disable_on_axis <- if (isTRUE(reindex_disable_on_axis)) { # TRUE
          bit64::seq.integer64(0L, ndim)
        } else if (isFALSE(reindex_disable_on_axis)) { # FALSE
          NULL
        } else { # NA
          ax <- bit64::seq.integer64(0L, ndim)
          ax[ax != self$axis]
        }
      }
      if (!is.null(reindex_disable_on_axis)) {
        stopifnot(
          "'reindex_disable_on_axis' must be a vector of integers" = (
            rlang::is_integerish(reindex_disable_on_axis) ||
              inherits(reindex_disable_on_axis, "integer64")
          ),
          "'reindex_disable_on_axis' must be finite" = is.finite(reindex_disable_on_axis),
          "'reindex_disable_on_axis' must be within the range of dimensions of the array" = all(
            reindex_disable_on_axis >= 0 & reindex_disable_on_axis <= ndim
          )
        )
        reindex_disable_on_axis <- unique(bit64::as.integer64(reindex_disable_on_axis))
      }
      private$.reindex_disable_on_axis <- reindex_disable_on_axis
      axes_to_reindex <- self$axes_to_reindex
      private$.reindexers <- vector("list", length = length(axes_to_reindex))
      shape <- self$array$shape()
      dnames <- self$array$dimnames()
      for (i in seq_along(axes_to_reindex)) {
        ax <- as.numeric(axes_to_reindex[i]) + 1L
        coords <- as.list(CoordsStrider$new(start = 0L, end = shape[ax] - 1L))
        coords <- if (length(coords) == 1L) {
          coords[[1L]]
        } else {
          unlist64(coords)
        }
        private$.reindexers[[i]] <- IntIndexer$new(coords)
        names(private$.reindexers)[i] <- dnames[ax]
      }
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
      private$.nextelems <- self$coords_axis$next_element()
      private$set_dim_points(dimnam, private$.nextelems)
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
    #' @field axes_to_reindex The axes to re-index
    #'
    axes_to_reindex = function() {
      ax <- bit64::seq.integer64(0L, self$array$ndim() - 1L)
      ax <- ax[!ax %in% self$reindex_disable_on_axis]
      if (length(ax)) {
        ax <- ax[ax != self$axis]
      }
      if (!length(ax)) {
        return(NULL)
      }
      return(ax)
    },
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
    reindex_disable_on_axis = function() private$.reindex_disable_on_axis,
    #' @field reindexable Shorthand to see if this iterator is poised to be
    #' re-indexed or not
    #'
    reindexable = function() length(self$axes_to_reindex) ||
      !bit64::as.integer64(self$axis) %in% self$reindex_disable_on_axis
  ),
  private = list(
    .array = NULL,
    .coords = list(),
    .axis = integer(1L),
    .nextelems = NULL,
    .reindex_disable_on_axis = NULL,
    .reindexers = list(),
    # @description Throw an error saying that re-indexed
    # iterators are not concatenatable
    .notConcatenatable = function() stop(errorCondition(
      message = "Re-indexed blockwise iterators are not concatenatable",
      class = "notConcatenatableError"
    )),
    # @description Reset internal state of SOMA Reader while keeping array open
    reset = function() {
      if (is.null(private$soma_reader_pointer)) {
        return(NULL)
      }
      sr_reset(private$soma_reader_pointer)
      return(invisible(NULL))
    },
    # @description Re-index an Arrow table
    reindex_arrow_table = function(tbl) {
      stopifnot(
        "'tbl' must be an Arrow table" = R6::is.R6(tbl) && inherits(tbl, 'Table')
      )
      dname <- self$array$dimnames()[self$axis + 1L]
      if (!dname %in% names(tbl)) {
        stop(
          "Cannot find ",
          sQuote(dname),
          " in the provided Arrow table",
          call. = FALSE
        )
      }
      op <- options(arrow.int64_downcast = FALSE)
      on.exit(options(op), add = TRUE, after = FALSE)
      coords <- self$coords
      coords[[dname]] <- CoordsStrider$new(
        private$.nextelems,
        stride = coords[[dname]]$stride
      )
      if (!bit64::as.integer64(self$axis) %in% self$reindex_disable_on_axis) {
        indexer <- IntIndexer$new(private$.nextelems)
        tbl[[dname]] <- indexer$get_indexer(
          tbl[[dname]]$as_vector(),
          nomatch_na = TRUE
        )
        rm(indexer)
      }
      for (dname in names(private$.reindexers)) {
        if (!dname %in% names(tbl)) {
          ""
        }
        indexer <- private$.reindexers[[dname]]
        tbl[[dname]] <- indexer$get_indexer(
          tbl[[dname]]$as_vector(),
          nomatch_na = TRUE
        )
      }
      attr(tbl, "coords") <- coords
      return(tbl)
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
    #' @description Concatenate the remainder of the blockwise iterator
    #'
    #' @return An Arrow Table with the remainder of the iterator
    #'
    concat = function() {
      if (self$reindexable) {
        private$.notConcatenatable()
      }
      return(soma_array_to_arrow_table_concat(self))
    }
  ),
  private = list(
    soma_reader_transform = function(x) {
      tbl <- soma_array_to_arrow_table(x)
      return(private$reindex_arrow_table(tbl))
    }
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
      reindex_disable_on_axis = NA
    ) {
      super$initialize(
        sr,
        array,
        coords,
        axis,
        ...,
        reindex_disable_on_axis = reindex_disable_on_axis
      )
      stopifnot(
        "Sparse reads only work with two-dimensional arrays" = self$array$ndim() == 2L
      )
      reprs <- c(
        'T',
        if (!bit64::as.integer64(0L) %in% self$reindex_disable_on_axis)'R',
        if (!bit64::as.integer64(1L) %in% self$reindex_disable_on_axis) 'C'
      )
      private$.repr <- match.arg(repr, choices = reprs)
      private$.shape <- sapply(coords, length)
    },
    #' @description Concatenate the remainder of the blockwise iterator
    #'
    #' @return A sparse matrix (determined by \code{self$repr}) with
    #' the remainder of the iterator
    #'
    concat = function() {
      if (self$reindexable) {
        private$.notConcatenatable()
      }
      return(soma_array_to_sparse_matrix_concat(self, private$.zero_based))
    }
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
    soma_reader_transform = function(x) {
      tbl <- private$reindex_arrow_table(soma_array_to_arrow_table(x))
      mat <- arrow_table_to_sparse(
        tbl,
        repr = self$repr,
        shape = private$.shape,
        zero_based = private$.zero_based
      )
      attr(mat, "coords") <- attr(tbl, "coords", exact = TRUE)
      return(mat)
    }
  )
)
