CoordsStrider <- R6::R6Class(
  classname = "CoordsStrider",
  cloneable = FALSE,
  public = list(
    initialize = function(coords, ..., stride = NULL, start = NULL, end = NULL) {
      if (missing(coords)) {
        stopifnot(
          rlang::is_integerish(start, 1L, TRUE) ||
            (inherits(start, "integer64") && length(start) == 1L && is.finite(start)),
          rlang::is_integerish(end, 1L, TRUE) ||
            (inherits(end, "integer64") && length(end) == 1L && is.finite(end)),
          start <= end
        )
        private$.start <- start
        private$.end <- end
        stride <- stride %||% abs(end - start + 1L)
        private$.index <- 0L
      } else {
        stopifnot(inherits(coords, c("integer64", "numeric", "integer")))
        private$.coords <- coords
        stride <- stride %||% length(coords)
        stopifnot(stride <= .Machine$integer.max)
        private$.index <- 1L
      }
      stopifnot(rlang::is_integerish(stride, 1L, TRUE) && stride > 0L)
      private$.stride <- stride
    },
    print = function() {
      cat("<", class(self)[1L], ">\n", sep = "")
      if (is.null(self$coords)) {
        cat("  start:", self$start, "\n")
        cat("  end:", self$end, "\n")
      } else {
        cat("  length(coords):", length(self$coords), "\n")
      }
      cat("  stride:", self$stride, "\n")
      return(invisible(self))
    },
    hasNext = function() {
      if (is.null(self$coords)) {
        return(private$.index <= abs(self$end - self$start))
      }
      return(private$.index <= length(self$coords))
    },
    nextElem = function() {
      if (!self$hasNext()) {
        private$.stopIteration()
      }
      if (is.null(self$coords)) {
        start <- min(self$start + private$.index, self$end)
        end <- min(start + self$stride - 1L, self$end)
        private$.index <- private$.index + self$stride
        if (start == end) {
          return(bit64::as.integer64(start))
        }
        by <- ifelse(start <= end, 1L, -1L)
        return(bit64::seq.integer64(from = start, to = end, by = by))
      }
      start <- private$.index
      end <- min(private$.index + self$stride - 1L, length(self$coords))
      private$.index <- end + 1
      return(self$coords[start:end])
    }
  ),
  active = list(
    coords = function() private$.coords,
    start = function() ifelse(is.null(self$coords), private$.start, min(self$coords)),
    end = function() ifelse(is.null(self$coords), private$.end, max(self$coords)),
    stride = function() private$.stride
  ),
  private = list(
    .coords = NULL,
    .start = NULL,
    .end = NULL,
    .stride = NULL,
    .index = NULL,
    .stopIteration = function() stop(errorCondition(
      "StopIteration",
      class = "stopIteration"
    ))
  )
)

#' @method as.list CoordsStrider
#' @export
#'
as.list.CoordsStrider <- function(x, ...) {
  f <- get('as.list.iter', envir = asNamespace('iterators'))
  return(f(x, ...))
}

#' @importFrom iterators nextElem
#' @export
#'
iterators::nextElem

#' @method nextElem CoordsStrider
#' @export
#'
nextElem.CoordsStrider <- function(obj, ...) {
  return(obj$nextElem())
}

#' @importFrom itertools hasNext
#' @export
#'
itertools::hasNext

#' @method hasNext CoordsStrider
#' @export
#'
hasNext.CoordsStrider <- function(obj, ...) {
  return(obj$hasNext())
}

unlist64 <- function(x) {
  stopifnot(
    is.list(x),
    all(vapply_lgl(x, inherits, what = 'integer64'))
  )
  res <- bit64::integer64(sum(vapply_int(x, length)))
  idx <- 1L
  for (i in seq_along(x)) {
    end <- idx + length(x[[i]])
    res[idx:(end - 1L)] <- x[[i]]
    idx <- end
  }
  return(res)
}
