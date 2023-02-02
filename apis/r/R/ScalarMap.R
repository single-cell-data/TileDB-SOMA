#' A Mapping Type for Scalars
#'
ScalarMap <- R6::R6Class(
  classname = 'ScalarMap',
  inherit = MappingBase,
  public = list(
    #' @param type ...
    #'
    #' @return ...
    initialize = function(type = 'any') {
      private$.type <- match.arg(arg = type, choices = .SCALAR_TYPES())
    },
    #' @param key ...
    #' @param value ...
    #'
    #' @return ...
    set = function(key, value) {
      stopifnot(is.null(x = value) || length(x = value) == 1L)
      if (private$.type != 'any') {
        stopifnot(inherits(x = value, what = private$.type))
      }
      super$set(key = key, value = value)
      return(invisible(x = self))
    }
  ),
  private = list(
    .type = character(length = 1L)
  )
)

.SCALAR_TYPES <- function() {
  return(c('any', 'numeric', 'integer', 'character', 'logical'))
}
