#' SOMA Read Iterator Base class
#'
#' Class that allows for read iteration of SOMA reads.
#' @keywords internal
#' @export

ReadIter <- R6::R6Class(
  classname = "ReadIter",

  public = list(

    #' @description Create (lifecycle: experimental)
    #' @param sr soma read pointer
    initialize = function(sr) {
      private$soma_reader_pointer <- sr
    },

    #' @description Check if iterated read is complete or not. (lifecycle: experimental)
    #' @return logical
    read_complete = function() {
      if (is.null(private$soma_reader_pointer)) {
          TRUE
      } else {
          sr_complete(private$soma_reader_pointer)
      }
    },

    #' @description Read the next chunk of an iterated read. (lifecycle: experimental).
    #' If read is complete, retunrs `NULL` and raises warning.
    #' @return \code{NULL} or one of arrow::\link[arrow]{Table}, \link{matrixZeroBasedView}
    read_next = function() {
      if (is.null(private$soma_reader_pointer)) {
          return(NULL)
      }
      if (self$read_complete()) {
        return(private$.readComplete())
      }
      rl <- sr_next(private$soma_reader_pointer)
      return(private$soma_reader_transform(rl))
    },

    #' @description  Concatenate remainder of iterator
    # to be refined in derived classes
    concat = function() {
      .NotYetImplemented()
    },

    #' @description Reset internal state of SOMA Reader while keeping array open
    reset = function() {
      if (is.null(private$soma_reader_pointer)) {
          return(NULL)
      }
      sr_reset(private$soma_reader_pointer)
    },

    #' @description Set dimension selection on given axis
    set_dim_point = function(dimname, points) {
      stopifnot("Name of dimension must be character" =
                    is.character(dimname),
                "Points must be int64 vector" =
                    is.vector(points) &&
                    inherits(points, "int64"))
      if (is.null(private$soma_reader_pointer)) {
          return(NULL)
      }
      sr_set_dim_points(private$soma_reader_pointer, dimname, points)
    }
  ),

  private = list(

    # Internal 'external pointer' object used for iterated reads
    soma_reader_pointer = NULL,

    # to be refined in derived classes
    soma_reader_transform = function(x) {
      .NotYetImplemented()
    },

    .readComplete = function(immediate. = FALSE) {
      warning(
        warningCondition(
          "Iteration complete, returning NULL",
          class = "iterationCompleteWarning"
        ),
        immediate. = immediate.
      )
      return(NULL)
    }

  )
)
