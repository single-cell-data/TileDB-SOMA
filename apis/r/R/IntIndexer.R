#' The SOMA Re-Indexer
#'
#' @description A re-indexer for unique integer indices
#' @export
IntIndexer <- R6::R6Class(
  classname = 'IntIndexer',
  public = list(
    #' @description Create a new re-indexer
    #'
    #' @param data Integer keys used to build the index (hash) table
    #'
    initialize = function(data) {
      stopifnot("'data' must be a vector of integers" = rlang::is_integerish(data, finite = TRUE) ||
                    (inherits(data, 'integer64') && all(is.finite(data))))

      # Setup the re-indexer with data
      private$.reindexer <- reindex_create()
      # Re-index
      reindex_map(private$.reindexer, bit64::as.integer64(data))
      return(invisible(NULL))
    },
    #' @description Get the underlying indices for the target data
    #'
    #' @param target Data to re-index
    #' @param nomatch The value to be returned when no match is found; will be
    #' coerced to a 64-bit integer
    #'
    #' @return A vector of 64-bit integers with \code{target} re-indexed
    #'
    get_indexer = function(target, nomatch = -1L) {
      # If `target` is an Arrow array, do Arrow handling
      if (R6::is.R6(target) && inherits(target, c('Array', 'ChunkedArray'))) {
        op <- options(arrow.int64_downcast = FALSE)
        on.exit(options(op), add = TRUE, after = FALSE)
        target <- target$as_vector()
        if (is.list(x = target)) {
          target <- unlist(x = target, use.names = FALSE)
        }
      }
      stopifnot(
        "'target' must be a vector or arrow array of integers" = rlang::is_integerish(target, finite = TRUE) ||
          (inherits(target, 'integer64') && all(is.finite(target))),
        "'nomatch' must be a single integer value" = rlang::is_integerish(x = nomatch, n = 1L)
      )
      # Do vector-based re-indexing
      val <- reindex_lookup(private$.reindexer, bit64::as.integer64(target))
      val[val == -1] <- bit64::as.integer64(nomatch)
      return(val)
    }
  ),
  private = list(
    # C++ reindexer
    .reindexer = NULL
  )
)
