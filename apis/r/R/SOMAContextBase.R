#' Base SOMA Context
#'
#' R6 mapping class for SOMA context options. This class should be used
#' as the basis for platform-specific contexts as it checks SOMA-specific
#' context options
#'
#' @keywords internal
#'
#' @export
#'
#' @noMd
#'
SOMAContextBase <- R6::R6Class(
  classname = 'SOMAContextBase',
  inherit = ScalarMap,
  public = list(
    #' @template param-config
    #'
    #' @return \Sexpr[results=rd]{tiledbsoma:::rd_return_virtual()}
    #'
    initialize = function(config = NULL) {
      calls <- vapply(
        X = lapply(X = sys.calls(), FUN = as.character),
        FUN = '[[',
        FUN.VALUE = character(length = 1L),
        1L
      )
      if ('SOMAContextBase$new' %in% calls) {
        stop(
            "'SOMAContextBase' is a virtual class and cannot be instantiated directly",
            call. = FALSE
        )
      }
      super$initialize()
      if (!is.null(x = config)) {
        msg <- "'config' must be a named vector"
        if (!is.vector(x = config)) {
          stop(msg)
        }
        conf_opts <- names(x = config)
        if (is.null(x = conf_opts) || !all(nzchar(x = conf_opts))) {
          stop(msg)
        }
        self$setv(config)
      }
    },
    #' @param key Key to set
    #' @templateVar key key
    #' @template param-value
    #'
    #' @return \[chainable\] Invisibly returns \code{self} with \code{value}
    #' added for \code{key}
    #'
    set = function(key, value) {
      super$set(key = key, value = value)
      soma_contexts <- .SOMA_CONTEXTS()
      if (key %in% names(x = soma_contexts)) {
        if (!inherits(x = private$.data[[key]], what = soma_contexts[key])) {
          stop(
            sQuote(x = key),
            " must be a ",
            sQuote(x = soma_contexts[key]),
            call. = FALSE
          )
        }
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
