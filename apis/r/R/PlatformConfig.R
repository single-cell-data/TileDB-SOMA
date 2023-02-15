#' Platform Configuration
#'
#' An R6 mapping type for configuring various \dQuote{options} for multiple
#' \dQuote{platforms}, essentially serves a multi-nested map where the inner
#' map is a \code{\link{ScalarMap}} contained within a \code{\link{ConfigList}}
#' (middle map):
#' \code{\{platform: \{op: \{key: value\}\}\}}
#'
#' @export
#'
PlatformConfig <- R6::R6Class(
  classname = 'PlatformConfig',
  inherit = MappingBase,
  public = list(
    #' @return The names of the \dQuote{platforms} (outer keys)
    #'
    platforms = function() {
      return(self$keys())
    },
    #' @param platform The \dQuote{platform} to pull option names (middle keys)
    #' for; pass \code{TRUE} to return all possible option names
    #'
    #' @return The option names (middle keys) for \code{platform}
    #'
    ops = function(platform = NULL) {
      platform <- platform %||% self$platforms()[1L]
      if (isTRUE(x = platform)) {
        ops <- Reduce(
          f = union,
          x = lapply(
            X = self$platforms(),
            FUN = \(p) private$.data[[p]]$keys()
          )
        )
        return(ops)
      }
      platform <- platform[1L]
      platform <- match.arg(arg = platform, choices = self$platforms())
      return(super$get(key = platform)$keys())
    },
    #' @param platform The name of the \dQuote{platform} (outer key) to fetch
    #' @param op The name of the \dQuote{option} of \code{platform} to fetch;
    #' if \code{NULL}, returns the \list[tiledbsoma:ConfigList]{configuration}
    #' for \code{platform}
    #' @param key The \dQuote{key} (inner key) for \code{op} in
    #' \code{platform} to fetch; if \code{NULL} and \code{op} is passed,
    #' returns the \link[tiledbsoma:ScalarMap]{map} for \code{op}
    #' in \code{platform}
    #' @templateVar key key
    #' @templateVar default null
    #' @template param-default
    #'
    #' @return The value of \code{key} for \code{op} in \code{platform} in the
    #' map, or \code{default} if \code{key} is not found
    #'
    get = function(platform, op = NULL, key = NULL, default = NULL) {
      platform <- platform[1L] %||% self$platforms()[1L]
      platform <- match.arg(arg = platform, choices = self$platforms())
      pmap <- super$get(key = platform)
      if (is.null(x = op)) {
        return(pmap)
      }
      return(pmap$get(op = op, key = key, default = default))
      # opmap <- pmap$get(key = op)
      # if (is.null(x = key)) {
      #   return(opmap)
      # }
      # return(opmap$get(key = key, default = default))
    },
    #' @param platform The name of the \dQuote{platform} (outer key) to fetch
    #'
    #' @return The \code{\link{ConfigList}} for \code{platform}
    #'
    get_ops = function(platform) {
      platform <- platform[1L] %||% self$platforms()[1L]
      platform <- match.arg(arg = platform, choices = self$platforms())
      return(super$get(key = platform))
    },
    #' @param platform The name of the \dQuote{platform} (outer key) to set
    #' @param op Name of the \dQuote{option} (middile key) in \code{platform}
    #' to set
    #' @param key Inner key to set
    #' @param value Value to add for \code{key}, or \code{NULL} to remove the
    #' entry for \code{key}; optionally provide only \code{platfomr}, \code{op},
    #' and \code{value} as a \code{\link{ScalarMap}} to update \code{op} for
    #' \code{platform} with the keys and values from \code{value}
    #'
    #' @return \[chainable\] Invisibly returns \code{self} with \code{value}
    #' added for \code{key} in \code{op}
    #'
    set = function(platform, op, key, value) {
      stopifnot(
        is.character(x = platform) && length(x = platform) == 1L,
        is.character(x = op) && length(x = op) == 1L
      )
      pmap <- super$get(key = platform, default = ConfigList$new())
      pmap$set(op = op, key = key, value = value)
      super$set(key = platform, value = pmap)
      return(invisible(x = self))
    },
    #' @template param-dots-ignored
    #'
    #' @return Nothing; \code{setv()} is disabled for
    #' \code{PlatformConfig} objects
    #'
    setv = function(...) {
      .NotYetImplemented()
    }
  )
)
