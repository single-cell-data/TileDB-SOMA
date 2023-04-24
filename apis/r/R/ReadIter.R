#' SOMA Read Iterator Base class
#'
#' Class that allows for read iteration of SOMA reads.

ReadIter <- R6::R6Class(
  classname = "ReadIter",

  public = list(
                
    uri = NULL,
    tiledbsoma_ctx = NULL,
    
    #' @description Create (lifecycle: experimental)
    initialize = function(uri, tiledbsoma_ctx) {
      self$uri <- uri
      self$tiledbsoma_ctx <- tiledbsoma_ctx
      private$soma_reader_setup()
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
          if (sr_complete(private$soma_reader_pointer)) {
              warning("Iteration complete, returning NULL")
              NULL
          } else {
              rl <- sr_next(private$soma_reader_pointer)
              private$soma_reader_transform(rl)
          }
      }
    },
    
    #' @description TODO
    # to be refined in derived classes
    concat = function() {
      NULL
    }

  ),

  private = list(

    # Internal 'external pointer' object used for iterated reads
    soma_reader_pointer = NULL,

    # Instantiate soma_reader_pointer with a soma_array_reader object
    soma_reader_setup = function() {
      private$soma_reader_pointer <- sr_setup(
        self$uri,
        config=as.character(tiledb::config(self$tiledbsoma_ctx$context()))
      )
    },

    ## to be refined in derived classes
    soma_reader_transform = function(x) {
      x
    }

  )
)
