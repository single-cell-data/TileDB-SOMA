#' @include ScalarList.R
#'
NULL

#' Base SOMA Context
#'
#' @keywords internal
#'
#' @export
#'
SOMAContextBase <- R6::R6Class(
  classname = 'SOMAContextBase',
  inherit = ScalarMap,
  public = list(
    #' @return ...
    initialize = function() {
      calls <- vapply(
        X = lapply(X = sys.calls(), FUN = as.character),
        FUN = '[[',
        FUN.VALUE = character(length = 1L),
        1L
      )
      if ('SOMAContextBase$new' %in% calls) {
        .NotYetImplemented()
      }
      super$initialize()
    },
    #' @param key ...
    #' @param value ...
    #'
    #' @return ...
    set = function(key, value) {
      super$set(key = key, value = value)
      soma_contexts <- .SOMA_CONTEXTS()
      if (key %in% names(x = soma_contexts)) {
        stopifnot(inherits(x = private$.data[[key]], what = soma_contexts[key]))
      }
      return(invisible(x = self))
    }
  )
)

#'
setClass(Class = 'SOMAContextBaseS4', contains = c('VIRTUAL', 'ScalarList'))

#' @method [ SOMAContextBaseS4
#' @export
#'
'[.SOMAContextBaseS4' <- function(x, i, ...) {
  .NotYetImplemented()
}

#' @method [[ SOMAContextBaseS4
#' @export
#'
"[[.SOMAContextBaseS4" <- function(x, i, ...) {
  i <- i[1L]
  i <- tryCatch(
    expr = match.arg(arg = i, choices = names(x = x)),
    error = \(...) i
  )
  soma_contexts <- .SOMA_CONTEXTS()
  val <- if (i %in% names(x = x)) {
    NextMethod()
  } else if (i %in% names(x = soma_contexts)) {
    vector(mode = soma_contexts[i])
  } else {
    NULL
  }
  if (length(x = val)) {
    names(x = val) <- i
  }
  return(val)
}

setMethod(
  f = '[<-',
  signature = c(x = 'SOMAContextBaseS4'),
  definition = function(x, i, ..., value) {
    .NotYetImplemented()
  }
)

setMethod(
  f = '[[<-',
  signature = c(x = 'SOMAContextBaseS4', i = 'ANY', j = 'missing', value = 'ANY'),
  definition = function(x, i, ..., value) {
    .NotYetImplemented()
  }
)

setMethod(
  f = '[[<-',
  signature = c(
    x = 'SOMAContextBaseS4',
    i = 'character',
    j = 'missing',
    value = 'ANY'
  ),
  definition = function(x, i, ..., value) {
    stopifnot(
      length(x = i) == length(x = value),
      length(x = i) == 1L,
      is.atomic(x = value)
    )
    slot(object = x, name = '.Data')[[i]] <- value
    return(x)
  }
)

setValidity(
  Class = 'SOMAContextBaseS4',
  method = function(object) {
    valid <- NULL
    browser()
    soma_contexts <- .SOMA_CONTEXTS()
    # Check individual entries
    for (i in names(x = object)) {
      if (i %in% names(x = soma_contexts)) {
        if (inherits(x = object[[i]], what = soma_contexts[[i]])) {
          valid <- c(
            valid,
            paste(
              "Context option",
              sQuote(x = i),
              "must be of type",
              sQuote(x = soma_contexts[[i]])
            )
          )
        }
      }
    }
    return(valid %||% TRUE)
  }
)

.SOMA_CONTEXTS <- function() {
  return(c(
    member_uris_are_relative = 'logical'
  ))
}
