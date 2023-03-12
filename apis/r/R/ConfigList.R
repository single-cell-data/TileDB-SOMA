#' A Configuration List
#'
#' An R6 mapping type for configuring various \dQuote{parameters}.
#' Essentially, serves as a nested map where the inner map is a
#' \code{\link{ScalarMap}}:
#' \code{\{<param>: \link[tiledbsoma:ScalarMap]{\{<key>: <value>\}}\}}
#'
#' @export
#'
#' @noMd
#'
ConfigList <- R6::R6Class(
  classname = 'ConfigList',
  inherit = MappingBase,
  public = list(
    #' @param param Outer key or \dQuote{parameter} to fetch
    #' @param key Inner key to fetch; pass \code{NULL} to return  the
    #' \link[tiledbsoma:ScalarMap]{map} for \code{param}
    #' @templateVar key key
    #' @templateVar default NULL
    #' @template param-default
    #'
    #' @return The value of \code{key} for \code{param} in the map, or
    #' \code{default} if \code{key} is not found
    #'
    get = function(param, key = NULL, default = quote(expr = )) {
      stopifnot(
        "'param' must be a single character value" = is_scalar_character(param)
      )
      parammap <- super$get(key = param, default = NULL)
      if (is.null(parammap)) {
        if (missing(default) || identical(x = default, y = quote(expr = ))) {
          private$.key_error(param)
        }
        return(default)
      }
      if (is.null(key)) {
        return(parammap)
      }
      return(parammap$get(key = key, default = default))
    },
    #' @param param Outer key or \dQuote{parameter} to set
    #' @param key Inner key to set
    #' @param value Value to add for \code{key}, or \code{NULL} to remove
    #' the entry for \code{key}; optionally provide only \code{param}
    #' and \code{value}
    #' as a \code{\link{ScalarMap}} to update \code{param} with the keys and
    #' values from \code{value}
    #'
    #' @return \[chainable\] Invisibly returns \code{self} with \code{value}
    #' added for \code{key} in \code{param}
    set = function(param, key, value) {
      stopifnot(
        "'param' must be a single character" = is_scalar_character(param)
      )
      parammap <- super$get(key = param, default = ScalarMap$new())
      if (missing(key) && inherits(x = value, what = 'ScalarMap')) {
        parammap$update(map = value)
        super$set(key = param, value = parammap)
        return(invisible(x = self))
      }
      parammap$set(key = key, value = value)
      super$set(key = param, value = parammap)
      return(invisible(self))
    },
    #' @template param-dots-ignored
    #'
    #' @return Nothing; \code{setv()} is disabled for \code{ConfigList} objects
    #'
    setv = function(...) {
      .NotYetImplemented()
    }
  )
)
