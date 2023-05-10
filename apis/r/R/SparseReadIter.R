#' SparseReadIter
#'
#' @description
#' `SparseReadIter` is a class that allows for iteration over 
#'  the results of a read operation from SOMA objects#' 
#' @export

SparseReadIter <- R6::R6Class(
  classname = "SparseReadIter",
  inherit = ReadIter,

  public = list(
                
    #' @description Create (lifecycle: experimental)
    initialize = function(uri, config, colnames = NULL, qc = NULL, dim_points = NULL, loglevel = "auto", repr) {
        # Initiate super class
          super$initialize (uri = uri, config = config, colnames = colnames, qc = qc,
                            dim_points = dim_points, loglevel = loglevel)
    
          private$repr <- repr
          
          # Get max soma dims for indeces
          tiledb_array <- tiledb::tiledb_array(uri)
          tiledb::tiledb_array_open(tiledb_array, type = "READ")
          max_soma_dim_0 <- as.integer(max(tiledb::tiledb_array_get_non_empty_domain_from_index(tiledb_array, 1)))
          max_soma_dim_1 <- as.integer(max(tiledb::tiledb_array_get_non_empty_domain_from_index(tiledb_array, 2)))
          tiledb::tiledb_array_close(tiledb_array)
          
          private$dims <- c(max_soma_dim_0, max_soma_dim_1)
    },
    
    ## refined from base class
    concat = function(){
  
      rl <- list()
      
      while (!self$read_complete()) {
        rl <- c(rl, self$read_next())
      }
      
      do.call(arrow::concat_tables, rl)
      
    }),
  
  private = list(
                 
    repr=NULL,
    dims=NULL,
    
    ## refined from base class
    soma_reader_transform = function(x) {
      arrow_table_to_sparse(soma_array_to_arrow(x), repr = private$repr, dims = private$dims)
    }
    
  )
)
