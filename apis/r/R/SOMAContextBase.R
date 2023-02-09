#' Base SOMA Context
#'
#' @keywords internal
#'
#' @export
#'
SOMAContextBase <- R6::R6Class(
  classname = 'SOMAContextBase',
  inherit = ScalarMap,
  public = list(
    #' @return ...
    initialize = function() {
      calls <- vapply(
        X = lapply(X = sys.calls(), FUN = as.character),
        FUN = '[[',
        FUN.VALUE = character(length = 1L),
        1L
      )
      if ('SOMAContextBase$new' %in% calls) {
        .NotYetImplemented()
      }
      super$initialize()
    },
    #' @param key ...
    #' @param value ...
    #'
    #' @return ...
    set = function(key, value) {
      super$set(key = key, value = value)
      soma_contexts <- .SOMA_CONTEXTS()
      if (key %in% names(x = soma_contexts)) {
        stopifnot(inherits(x = private$.data[[key]], what = soma_contexts[key]))
      }
      return(invisible(x = self))
    }
  )
)

.SOMA_CONTEXTS <- function() {
  return(c(
    member_uris_are_relative = 'logical'
  ))
}
