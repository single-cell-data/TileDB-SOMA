#' The Coordinate Strider
#'
#' @description The \code{CoordsStrider} allows creating coordinate slices
#' in an interated manner. Alternatively, it can chunk an existing vector of
#' coordinates
#'
#' @note The \code{CoordsStrider} operates using
#' \link[bit64:integer64]{64-bit integer} objects; as such, accessing fields,
#' such as \code{strider$start} or \code{strider$stride} will return an
#' \code{integer64} object, which functions differently than a regular
#' \code{integer}. Use with caution and convert back to integers or numerics
#' as necessary
#'
#' @keywords internal
#'
#' @export
#'
#' @examples
#' strider <- CoordsStrider$new(start = 1L, end = 200L, stride = 60L)
#' while (strider$has_next()) {
#'   str(strider$next_element())
#' }
#'
CoordsStrider <- R6::R6Class(
  classname = "CoordsStrider",
  cloneable = FALSE,
  public = list(
    #' @description Create a coordinate strider
    #'
    #' @param coords An integer vector of coordinates
    #' @param ... Ignored
    #' @param stride The stride of how many coordinates to yield per iteration;
    #' by default, will try to yield all coordinates per iteration
    #' @param start If \code{coords} is missing, the starting coordinate
    #' to generate
    #' @param end If \code{coords} is missing, the ending coordinate
    #' to generate
    #'
    initialize = function(coords, ..., stride = NULL, start = NULL, end = NULL) {
      if (missing(coords)) {
        stopifnot(
          "'start' must be a single integer value" = rlang::is_integerish(start, 1L, TRUE) ||
            (inherits(start, "integer64") && length(start) == 1L && is.finite(start)),
          "'end' must be a single integer value"  = rlang::is_integerish(end, 1L, TRUE) ||
            (inherits(end, "integer64") && length(end) == 1L && is.finite(end)),
          "'start' must be less than or equal to 'end'" = start <= end
        )
        private$.start <- bit64::as.integer64(start)
        private$.end <- bit64::as.integer64(end)
        stride <- stride %||% bit64::abs.integer64(self$end - self$start + 1L)
        private$.index <- 0L
      } else {
        stopifnot(
          "'coords' must be a vector of integer-like values" = inherits(coords, c("integer64", "numeric", "integer"))
        )
        private$.coords <- bit64::as.integer64(coords)
        stride <- stride %||% length(coords)
        stopifnot(
          "'stride' must be less than `Machine$integer.max`" = stride <= .Machine$integer.max
        )
        private$.index <- 1L
      }
      stopifnot(
        "'stride' must be a single integer value" = rlang::is_integerish(stride, 1L, TRUE) ||
          (inherits(stride, "integer64") && length(stride == 1L) && is.finite(stride)),
        stride > 0L
      )
      private$.stride <- bit64::as.integer64(stride)
    },
    #' @description Print the coordinate strider to the screen
    #'
    print = function() {
      cat("<", class(self)[1L], ">\n", sep = "")
      if (is.null(self$coords)) {
        cat("  start:", format(self$start), "\n")
        cat("  end:", format(self$end), "\n")
      } else {
        cat("  length(coords):", length(self$coords), "\n")
      }
      cat("  stride:", format(self$stride), "\n")
      return(invisible(self))
    },
    #' @description Get the length (span) of the coordinates within the strider
    #'
    #' @return The length (span) of the coordinate strider
    #'
    length = function() {
      if (is.null(self$coords)) {
        len <- as.numeric(abs(self$end - self$start))
        if (!self$start) {
          len <- len + 1L
        }
        return(len)
      }
      return(length(self$coords))
    },
    #' @description Determine if there are more coordinates to yield
    #'
    #' @return \code{TRUE} if there are more coordinates to yield or
    #' \code{FALSE} if otherwise
    #'
    has_next = function() {
      if (is.null(self$coords)) {
        return(private$.index <= abs(self$end - self$start))
      }
      return(private$.index <= length(self$coords))
    },
    #' @description Generate the next set of coordinates to yield. If there are
    #' no more coordinates to yield, raises a \code{stopIteration} error
    #'
    #' @return An integer vector of the next set of coordinates
    #'
    next_element = function() {
      if (!self$has_next()) {
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
      end <- as.integer(bit64::min.integer64(
        private$.index + self$stride - 1L,
        length(self$coords)
      ))
      private$.index <- end + 1L
      return(self$coords[start:end])
    }
  ),
  active = list(
    #' @field coords If set, the coordinates to iterate over
    #'
    coords = function() private$.coords,
    #' @field start If set, the starting point of the iterated coordinates;
    #' otherwise the minimum value of \code{self$coords}
    #'
    start = function() if (is.null(self$coords)) {
      private$.start
    } else {
      bit64::min.integer64(self$coords)
    },
    #' @field end If set, the end point of the iterated coordinates;
    #' otherwise the maximum value of \code{self$coords}
    #'
    end = function() if (is.null(self$coords)) {
      private$.end
    } else {
      bit64::max.integer64(self$coords)
    },
    #' @field stride The stride, or how many coordinates to generate per
    #' iteration; note: this field is settable, which will reset the iterator
    #'
    stride = function(value) {
      if (missing(value)) {
        return(private$.stride)
      }
      stopifnot(
        "'stride' must be a single integer value" = (rlang::is_integerish(value, n = 1L, finite = TRUE) ||
          (inherits(value, "integer64") && length(value) == 1L && is.finite(value))
        ) &&
          value > 0L
      )
      private$.stride <- bit64::as.integer64(value)
      index <- ifelse(is.null(self$coords), yes = 0L, no = 1L)
      if (private$.index != index) {
        warning(warningCondition(
          "The stride has changed, coordinates have been put back at the beginning",
          class = "coordsStrideChangedWarning"
        ))
      }
      private$.index <- index
      return(invisible(NULL))
    }
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
  res <- vector(mode = "list", length = ceiling(x$length() / x$stride))
  i <- 1L
  while (x$has_next()) {
    if (i > length(res)) {
      stop("this shouldn't be reached")
    }
    res[i] <- list(x$next_element())
    i <- i + 1L
  }
  return(Filter(Negate(is.null), res))
}

#' @method length CoordsStrider
#' @export
#'
length.CoordsStrider <- function(x) x$length()

#' @exportS3Method iterators::nextElem
#'
nextElem.CoordsStrider <- function(obj, ...) obj$next_element()

#' @exportS3Method itertools::hasNext
#'
hasNext.CoordsStrider <- function(obj, ...) obj$has_next()

unlist64 <- function(x) {
  stopifnot(
    "'x' must be a list" = is.list(x),
    "'x' must contain 'integer64' values" = all(vapply_lgl(x, inherits, what = 'integer64'))
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
