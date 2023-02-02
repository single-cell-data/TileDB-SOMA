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
    #'
    initialize = function(...) {
      calls <- vapply(
        X = lapply(X = sys.calls(), FUN = as.character),
        FUN = '[[',
        FUN.VALUE = character(length = 1L),
        1L
      )
      if ('MappingBase$new' %in% calls) {
        .NotYetImplemented()
      }
      private$.data <- list()
    },
    #' @return The keys of the map
    #'
    keys = function() {
      return(names(x = private$.data))
    },
    #' @return The values of the map
    #'
    values = function() {
      return(unname(obj = self$items()))
    },
    #' @return Return the items of the map as a list
    #'
    items = function() {
      return(private$.data)
    },
    #' @param key Key to get
    #' @param default Default value to return if \code{key} is not found
    #'
    #' @return The value of \code{key} in the map, or \code{default} if
    #' \code{key} is not found
    #'
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
    #' @param key Key to set
    #' @param value Value to add for \code{key}
    #'
    #' @return \[chainable\] Invisibly returns \code{self} with
    #' \code{value} added as \code{key}
    #'
    set = function(key, value) {
      stopifnot(
        is.character(x = key),
        # is.null(x = value) || length(x = key) == length(x = value),
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
    #' @param ... Named arguments to add to \code{self}
    #'
    #' @return \[chainable\] Invisibly returns \code{self} with the values
    #' of \code{...}
    #'
    setv = function(...) {
      args <- as.list(x = c(...))
      stopifnot(is_named_list(x = args))
      for (i in seq_along(along.with = args)) {
        self$set(key = names(x = args)[i], value = args[[i]])
      }
      return(invisible(x = self))
    },
    #' @param map A mapping type to update the current map with
    #' @return \[chainable\] Invisibly returns \code{self} with the value
    #' of \code{map}
    update = function(map) {
      stopifnot(inherits(x = map, what = 'MappingBase'))
      self$setv(map$items())
      return(invisible(x = self))
    },
    #' @return The number of items in the map
    length = function() {
      return(length(x = private$.data))
    },
    #' @return The map as a list
    #'
    to_list = function() {
      return(private$.data)
    },
    #' @return \[chainable\] Prints information about the map to \code{stdout}
    #' and invisibly returns \code{self}
    print = function() {
      vowels <- c('a', 'e', 'i', 'o', 'u')
      cls <- class(x = self)[1L]
      prd <- ifelse(
        test = tolower(x = substr(x = cls, start = 1L, stop = 2L)) %in% vowels,
        yes = 'An',
        no = 'A'
      )
      lng <- self$length()
      ent <- ifelse(test = lng == 1L, yes = 'entry', no = 'entries')
      cat(prd, cls, "map with", lng, ent, sep = ' ')
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
