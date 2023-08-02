#' Base SOMA Context
#'
#' R6 mapping class for SOMA context options. This class should be used
#' as the basis for platform-specific contexts as it checks SOMA-specific
#' context options
#'
#' @keywords internal
#' @export

SOMAContextBase <- R6::R6Class(
  classname = 'SOMAContextBase',
  inherit = ScalarMap,
  public = list(
    #' @template param-config
    #'
    #' @return \Sexpr[results=rd]{tiledbsoma:::rd_return_virtual()}
    #'
    initialize = function(config = NULL) {
      calls <- vapply_char(
        X = lapply(X = sys.calls(), FUN = as.character),
        FUN = '[[',
        1L
      )
      if ('SOMAContextBase$new' %in% calls) {
        stop(
            "'SOMAContextBase' is a virtual class and cannot be instantiated directly",
            call. = FALSE
        )
      }
      super$initialize()
      if (!is.null(config)) {
        msg <- "'config' must be a named vector"
        if (!is.vector(config)) {
          stop(msg, call. = FALSE)
        }
        if (!is_named(config, allow_empty = FALSE)) {
          stop(msg, call. = FALSE)
        }
        # conf_opts <- names(config)
        # if (is.null(conf_opts) || !all(nzchar(conf_opts))) {
        #   stop(msg)
        # }
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
      if (key %in% names(soma_contexts)) {
        if (!inherits(x = private$.data[[key]], what = soma_contexts[key])) {
          stop(
            sQuote(key),
            " must be a ",
            sQuote(soma_contexts[key]),
            call. = FALSE
          )
        }
      }
      return(invisible(self))
    }
  )
)

.SOMA_CONTEXTS <- function() {
  return(c(
    member_uris_are_relative = 'logical'
  ))
}
