#' SOMA Read Iterator Base class
#'
#' Class that allows for read iteration of SOMA reads.

ReadIter <- R6::R6Class(
  classname = "ReadIter",

  public = list(
                
    #' @description Create (lifecycle: experimental)
    initialize = function(uri, config, colnames = NULL, qc = NULL, dim_points = NULL, loglevel = "auto") {
        # Instantiate soma_reader_pointer with a soma_array_reader object
          private$soma_reader_pointer <- sr_setup(
            uri = uri,
            config = config,
            colnames = colnames,
            qc = qc,
            dim_points = dim_points,
            loglevel = loglevel
          )
    },

    #' @description Check if iterated read is complete or not. (lifecycle: experimental)
    read_complete = function() {
      if (is.null(private$soma_reader_pointer)) {
          TRUE
      } else {
          sr_complete(private$soma_reader_pointer)
      }
    },

    #' @description Read the next chunk of an iterated read. (lifecycle: experimental)
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
    
    #' @description TODO
    # to be refined in derived classes
    concat = function() {
      .NotYetImplemented()
    }

  ),

  private = list(

    # Internal 'external pointer' object used for iterated reads
    soma_reader_pointer = NULL,

    ## to be refined in derived classes
    soma_reader_transform = function(x) {
      .NotYetImplemented()
    }

  )
)
