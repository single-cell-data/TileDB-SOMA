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
      private$.reindexer <- CPP_INDEXER(private$.context$context())
      CPP_MAP_LOCATIONS(private$.reindexer, data)
      return(invisible(NULL))
    },
    #' @description ...
    #'
    #' @param target Data to re-index
    #'
    #' @return ...
    #'
    get_indexer = function(target) {
      stopifnot(
        "'target' must be a vector or arrow array of integers" = rlang::is_integerish(
          target,
          finite = TRUE
        ) ||
          (inherits(target, 'integer64') && all(is.finite(data))) ||
          (R6::is.R6(target) && inherits(target, c('Array', 'ChunkedArray')))
      )
      if (R6::is.R6(target) && inherits(target, c('Array', 'ChunkedArray'))) {
        return(CPP_GET_INDEXER_ARROW(private$.reindexer, target))
      }
      return(CPP_GET_INDEXER_GENERAL(private$.reindexer, target))
    }
  ),
  private = list(
    .context = NULL,
    .reindexer = NULL
  )
)
