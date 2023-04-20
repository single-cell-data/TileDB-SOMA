#' Platform Configuration
#'
#' An R6 mapping type for configuring various \dQuote{parameters} for multiple
#' \dQuote{platforms}, essentially serves a multi-nested map where the inner
#' map is a \code{\link{ScalarMap}} contained within a \code{\link{ConfigList}}
#' (middle map):
#' \code{\{platform: \{param: \{key: value\}\}\}}
#'
#' @export
#'
#' @noMd
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
    #' @param platform The \dQuote{platform} to pull parameter names
    #' (middle keys)
    #' for; pass \code{TRUE} to return all possible parameter names
    #'
    #' @return The parameter names (middle keys) for \code{platform}
    #'
    params = function(platform = NULL) {
      stopifnot(
        "'platform' must be a scalar character or logical value" = is.null(platform) ||
          is_scalar_character(platform) ||
          is_scalar_logical(platform)
      )
      platform <- platform %||% self$platforms()[1L]
      if (isTRUE(x = platform)) {
        params <- Reduce(
          f = union,
          x = lapply(
            X = self$platforms(),
            FUN = function(p) {
              private$.data[[p]]$keys()
            }
          )
        )
        return(params)
      }
      platform <- match.arg(arg = platform, choices = self$platforms())
      return(super$get(key = platform)$keys())
    },
    #' @param platform The name of the \dQuote{platform} (outer key) to fetch
    #' @param param The name of the \dQuote{paramters} of \code{platform}
    #' to fetch; if \code{NULL}, returns the
    #' \link[tiledbsoma:ConfigList]{configuration} for \code{platform}
    #' @param key The \dQuote{key} (inner key) for \code{param} in
    #' \code{platform} to fetch; if \code{NULL} and \code{param} is passed,
    #' returns the \link[tiledbsoma:ScalarMap]{map} for \code{param}
    #' in \code{platform}
    #' @templateVar key key
    #' @templateVar default null
    #' @template param-default
    #'
    #' @return The value of \code{key} for \code{param} in \code{platform} in the
    #' map, or \code{default} if \code{key} is not found
    #'
    get = function(
      platform,
      param = NULL,
      key = NULL,
      default = quote(expr = )
    ) {
      if (!length(self)) {
        warning("No platforms configured", call. = FALSE)
        return(NULL)
      }
      stopifnot(
        "'platform' must be a single character value" = is.null(platform) ||
          is_scalar_character(platform)
      )
      platform <- platform %||% self$platforms()[1L]
      pmap <- super$get(key = platform, default = NULL)
      if (is.null(pmap)) {
        if (missing(default) || identical(x = default, y = quote(expr = ))) {
          private$.key_error(platform)
        }
        return(default)
      }
      if (is.null(param)) {
        return(pmap)
      }
      return(pmap$get(param = param, key = key, default = default))
    },
    #' @param platform The name of the \dQuote{platform} (outer key) to fetch
    #'
    #' @return The \code{\link{ConfigList}} for \code{platform}
    #'
    get_params = function(platform) {
      stopifnot(
        "'platform' must be a single character value" = is.null(platform) ||
          is_scalar_character(platform)
      )
      platform <- platform %||% self$platforms()[1L]
      platform <- match.arg(arg = platform, choices = self$platforms())
      return(super$get(platform))
    },
    #' @param platform The name of the \dQuote{platform} (outer key) to set
    #' @param param Name of the \dQuote{parameter} (middle key) in
    #' \code{platform} to set
    #' @param key Inner key to set
    #' @param value Value to add for \code{key}, or \code{NULL} to remove the
    #' entry for \code{key}; optionally provide only \code{platfomr},
    #' \code{param}, and \code{value} as a \code{\link{ScalarMap}} to
    #' update  \code{param} for \code{platform} with the keys and values
    #' from \code{value}
    #'
    #' @return \[chainable\] Invisibly returns \code{self} with \code{value}
    #' added for \code{key} in \code{param} for \code{platform}
    #'
    set = function(platform, param, key, value) {
      stopifnot(
        "'platform' must be a single character value" = is_scalar_character(platform),
        "'param' must be a single character value" = is_scalar_character(param)
      )
      pmap <- super$get(key = platform, default = ConfigList$new())
      pmap$set(param = param, key = key, value = value)
      super$set(key = platform, value = pmap)
      return(invisible(self))
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
