#' SOMA Read Iterator Base class
#'
#' Dispatcher class to read full SOMAArray

ReadFull <- R6::R6Class(
  classname = "ReadFull",

  public = list(
                
    #' @description Create (lifecycle: experimental)
    initialize = function(uri, config, colnames = NULL, qc = NULL, dim_points = NULL, loglevel = "auto") {
        # Instantiate soma_reader_pointer with a soma_array_reader object
          private$rl <- soma_array_reader(
            uri = uri,
            config = config,
            colnames = colnames,
            qc = qc,
            dim_points = dim_points,
            loglevel = loglevel
          )
    },

    #' @description Check if iterated read is complete or not. (lifecycle: experimental)
    read = function() {
      return(private$soma_reader_transform(private$rl))
    }

  ),

  private = list(

    # Internal 'external pointer' object used for iterated reads
    rl = NULL,

    ## to be refined in derived classes
    soma_reader_transform = function(x) {
      .NotYetImplemented()
    }

  )
)
