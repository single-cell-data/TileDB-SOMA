#' SOMA Read Iterator Base class
#'
#' Virtual class that allows for read iteration of SOMA reads.
#'
#' @keywords internal
#'
#' @export
#'
#' @seealso \code{\link{BlockwiseReadIterBase}},
#' \code{\link{SparseReadIter}},
#' \code{\link{TableReadIter}}
#'
ReadIter <- R6::R6Class(
  classname = "ReadIter",
  public = list(

    #' @description Create (lifecycle: maturing).
    #'
    #' @param sr soma read pointer.
    #'
    initialize = function(sr) {
      private$soma_reader_pointer <- sr
    },

    #' @description @description Check if iterated read is complete or not
    #' (lifecycle: maturing).
    #'
    #' @return logical
    #'
    read_complete = function() {
      if (is.null(private$soma_reader_pointer)) {
        TRUE
      } else {
        mq_complete(private$soma_reader_pointer)
      }
    },

    #' @description Read the next chunk of an iterated read. If read is
    #' complete, returns \code{NULL} and raises warning (lifecycle: maturing).
    #'
    #' @return \code{NULL} or one of \link[arrow:Table]{arrow::Table},
    #' \link{matrixZeroBasedView}.
    #'
    read_next = function() {
      if (is.null(private$soma_reader_pointer)) {
        return(NULL)
      }
      if (self$read_complete()) {
        return(private$.readComplete())
      }

      return(private$.read_next())
    },

    #' @description  Concatenate remainder of iterator.
    # to be refined in derived classes
    #
    concat = function() {
      .NotYetImplemented()
    }
  ),
  private = list(

    # Internal 'external pointer' object used for iterated reads
    #
    soma_reader_pointer = NULL,

    # to be refined in derived classes
    #
    soma_reader_transform = function(x) .NotYetImplemented(),

    # Internal `read_next()` to avoid `self$read_complete()` checks
    #
    .read_next = function() {
      if (is.null(private$soma_reader_pointer)) {
        return(NULL)
      }
      rl <- mq_next(private$soma_reader_pointer)
      return(private$soma_reader_transform(rl))
    },

    # Throw a warning for read completion
    #
    .readComplete = function() {
      warning(warningCondition(
        "Iteration complete, returning NULL",
        class = "iterationCompleteWarning"
      ))
      return(NULL)
    }
  )
)
