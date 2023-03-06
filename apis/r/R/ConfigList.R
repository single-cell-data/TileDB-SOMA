#' A Configuration List
#'
#' An R6 mapping type for configuring various \dQuote{options}. Essentially,
#' serves as a nested map where the inner map is a \code{\link{ScalarMap}}:
#' \code{\{<op>: \link[tiledbsoma:ScalarMap]{\{<key>: <value>\}}\}}
#'
#' @export
#'
#' @noMd
#'
ConfigList <- R6::R6Class(
  classname = 'ConfigList',
  inherit = MappingBase,
  public = list(
    #' @param op Outer key or \dQuote{option} to fetch
    #' @param key Inner key to fetch; pass \code{NULL} to return  the
    #' \link[tiledbsoma:ScalarMap]{map} for \code{op}
    #' @templateVar key key
    #' @templateVar default NULL
    #' @template param-default
    #'
    #' @return The value of \code{key} for \code{op} in the map, or
    #' \code{default} if \code{key} is not found
    #'
    get = function(op, key = NULL, default = rlang::missing_arg()) {
      op <- op[1L]
      op <- tryCatch(
        expr = match.arg(arg = op, choices = self$keys()),
        error = \(...) NULL
      )
      if (is.null(x = op)) {
        if (rlang::is_missing(x = default)) {
          private$.key_error(key = op)
        }
        return(default)
      }
      opmap <- super$get(key = op)
      if (is.null(x = key)) {
        return(opmap)
      }
      return(opmap$get(key = key, default = default))
    },
    #' @param op Outer key or \dQuote{option} to set
    #' @param key Inner key to set
    #' @param value Value to add for \code{key}, or \code{NULL} to remove the
    #' entry for \code{key}; optionally provide only \code{op} and \code{value}
    #' as a \code{\link{ScalarMap}} to update \code{op} with the keys and
    #' values from \code{value}
    #'
    #' @return \[chainable\] Invisibly returns \code{self} with \code{value}
    #' added for \code{key} in \code{op}
    set = function(op, key, value) {
      stopifnot(
        "'op' must be a single character" = is_scalar_character(op)
      )
      opmap <- super$get(key = op, default = ScalarMap$new())
      if (missing(x = key) && inherits(x = value, what = 'ScalarMap')) {
        opmap$update(map = value)
        super$set(key = op, value = opmap)
        return(invisible(x = self))
      }
      stopifnot(is_scalar_character(op))
      opmap$set(key = key, value = value)
      super$set(key = op, value = opmap)
      return(invisible(x = self))
    },
    #' @template param-dots-ignored
    #'
    #' @return Nothing; \code{setv()} is disabled for \code{ConfigList} objects
    #'
    setv = function(...) {
      .NotYetImplemented()
    }
    # add_map = function(key, value) {
    #   stopifnot(
    #     is.character(x = key),
    #     length(x = key) == 1L,
    #     inherits(x = value, what = 'ScalarMap')
    #   )
    #   private$.data[[key]] <- value
    #   return(invisible(x = self))
    # },
    # set_from_list = function(key, value) {
    #   stopifnot(is_named_list(x = value))
    #   cfg <- ScalarMap$new()
    #   for (i in names(x = value)) {
    #     cfg$set(key = i, value = value[[i]])
    #   }
    #   self$add_map(key = key, value = cfg)
    #   return(invisible(x = self))
    # }
  )
)
