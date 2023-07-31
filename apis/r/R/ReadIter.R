#' SOMA Read Iterator Base class
#'
#' Class that allows for read iteration of SOMA reads.

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
          NULL
      } else {
          if (self$read_complete()) {
              warning("Iteration complete, returning NULL")
              NULL
          } else {
              rl <- sr_next(private$soma_reader_pointer)
              return(private$soma_reader_transform(rl))
          }
      }
    },
    
    #' @description  Concatenate remainder of iterator
    # to be refined in derived classes
    concat = function() {
      .NotYetImplemented()
    }

  ),

  private = list(

    # Internal 'external pointer' object used for iterated reads
    soma_reader_pointer = NULL,

    # to be refined in derived classes
    soma_reader_transform = function(x) {
      .NotYetImplemented()
    }

  )
)
