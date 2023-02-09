#' Platform Configuration
#'
#' @export
#'
PlatformConfig <- R6::R6Class(
  classname = 'PlatformConfig',
  inherit = MappingBase,
  public = list(
    #' @return ...
    platforms = function() {
      return(self$keys())
    },
    #' @param platform ...
    #'
    #' @return ...
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
    #' @param platform ...
    #' @param op ...
    #' @param key ...
    #' @param default ...
    #'
    #' @return ...
    get = function(platform, op = NULL, key = NULL, default = NULL) {
      platform <- platform[1L] %||% self$platforms()[1L]
      platform <- match.arg(arg = platform, choices = self$platforms())
      pmap <- super$get(key = platform)
      if (is.null(x = op)) {
        return(pmap)
      }
      opmap <- pmap$get(key = op)
      if (is.null(x = key)) {
        return(opmap)
      }
      return(opmap$get(key = key, default = default))
    },
    #' @param platform ...
    #'
    #' @return ...
    get_ops = function(platform) {
      platform <- platform[1L] %||% self$platforms()[1L]
      platform <- match.arg(arg = platform, choices = self$platforms())
      return(super$get(key = platform))
    },
    #' @param platform ...
    #' @param op ...
    #' @param key ...
    #' @param value ...
    #'
    #' @return ...
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
    setv = function(...) {
      .NotYetImplemented()
    }
  )
)
