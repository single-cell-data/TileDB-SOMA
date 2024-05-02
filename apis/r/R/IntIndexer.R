#' The SOMA Re-Indexer
#'
#' @description A re-indexer for unique integer indices
#'
IntIndexer <- R6::R6Class(
  classname = 'IntIndexer',
  public = list(
    #' @description Create a new re-indexer
    #' @param data Integer keys used to build the index (hash) table
    #' @param tiledbsoma_ctx Optional \code{\link{SOMATileDBContext}}
    #'
    initialize = function(data, tiledbsoma_ctx = NULL) {
      stopifnot(
        "'data' must be a vector of integers" = rlang::is_integerish(data, finite = TRUE) ||
          (inherits(data, 'integer64') && all(is.finite(data))),
        "'tiledbsoma_ctx' must be a SOMATileDBContext" = is.null(tiledbsoma_ctx) ||
          inherits(tiledbsoma_ctx, 'SOMATileDBContext')
      )
      private$.context <- tiledbsoma_ctx %||% SOMATileDBContext$new()
      # Setup the re-indexer
      private$.reindexer <- CPP_INDEXER(private$.context$context())
      # Setup the keys for the reindexer
      CPP_MAP_LOCATIONS(private$.reindexer, bit64::as.integer64(data))
      return(invisible(NULL))
    },
    #' @description Get the underlying indices for the target data
    #'
    #' @param target Data to re-index
    #'
    #' @return A vector of 64-bit integers with \code{target} re-indexed
    #'
    get_indexer = function(target) {
      stopifnot(
        "'target' must be a vector or arrow array of integers" = rlang::is_integerish(target, finite = TRUE) ||
          (inherits(target, 'integer64') && all(is.finite(data))) ||
          (R6::is.R6(target) && inherits(target, c('Array', 'ChunkedArray')))
      )
      # If `target` is an Arrow array, do Arrow handling
      if (R6::is.R6(target) && inherits(target, c('Array', 'ChunkedArray'))) {
        return(CPP_GET_INDEXER_ARROW(private$.reindexer, target))
      }
      # Do vector-based re-indexing
      # Cast to integer64 for C++
      return(bit64::as.integer64(CPP_GET_INDEXER_GENERAL(
        private$.reindexer,
        bit64::as.integer64(target)
      )))
    }
  ),
  private = list(
    # TileDB SOMA Context
    .context = NULL,
    # C++ reindexer
    .reindexer = NULL
  )
)
