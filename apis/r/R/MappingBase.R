#' R6 Base Mapping Type
#'
#' @keywords internal
#'
#' @export
#'
MappingBase <- R6::R6Class(
  classname = 'MappingBase',
  lock_class = TRUE,
  public = list(
    #' @param ... Ignored
    #'
    #' @return ...
    initialize = function(...) {
      .NotYetImplemented()
    },
    #' @return The keys of the map
    keys = function() {
      return(names(x = private$.data))
    },
    #' @return The values of the map
    values = function() {
      return(unname(obj = self$items()))
    },
    #' @return ...
    items = function() {
      return(private$.data)
    },
    #' @param key ...
    #' @param default ...
    #'
    #' @return ...
    get = function(key, default = NULL) {
      key <- tryCatch(
        expr = match.arg(arg = key, choices = self$keys()),
        error = \(...) NULL
      )
      if (is.null(x = key)) {
        return(default)
      }
      return(private$.data[[key]])
    },
    #' @param key ...
    #' @param value ...
    #'
    #' @return \[chainable\] Invisibly returns \code{self}
    set = function(key, value) {
      stopifnot(
        is.character(x = key),
        !is.null(x = value) && length(x = key) == length(x = value),
        length(x = key) == 1L
      )
      private$.data[[key]] <- value
      # TODO: figure out how to get x$`key` to act as active binding
      # if (is.null(x = value)) {
      #   self$.__enclos_env__$.__active__[[key]] <- NULL
      # } else {
      #   f <- paste(
      #     'function(value) {',
      #     '  if (missing(x = value)) {',
      #     paste0('    return(private$.data[[', sQuote(x = key, q = FALSE), ']])'),
      #     '  }',
      #     paste0('  private$.data[[', sQuote(x = key, q = FALSE), ']] <- value'),
      #     '  return(invisible(x = NULL))',
      #     '}',
      #     sep = '\n'
      #   )
      #   f <- eval(expr = str2expression(text = f), envir = self$.__enclos_env__)
      #   makeActiveBinding(sym = key, fun = f, env = self$.__enclos_env__)
      #   # self$.__enclos_env__$.__active__[[key]] <- f
      # }
      private$.data <- Filter(f = length, x = private$.data)
      # self <- self$clone(deep = TRUE)
      return(invisible(x = self))
    },
    #' @param ...
    #'
    #' @return \[chainable\] Invisibly returns \code{self}
    setv = function(...) {
      args <- list(...)
      stopifnot(is_named_list(x = args))
      for (i in seq_along(along.with = args)) {
        self$set(key = names(x = args)[i], value = args[[i]])
      }
      return(invisible(x = self))
    }
  ),
  private = list(
    .data = list()
  )
)

#' @method [[ MappingBase
#' @export
#'
'[[.MappingBase' <- function(x, i, ..., default = NULL) {
  return(x$get(key = i, default = default))
}


#' @method [[<- MappingBase
#' @export
#'
'[[<-.MappingBase' <- function(x, i, ..., value) {
  stopifnot(is.character(x = i), length(x = i) == 1L)
  x$set(key = i, value = value)
  return(x)
}
