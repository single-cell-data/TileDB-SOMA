#' @include MappingBase.R
#'
NULL

ConfigList <- R6::R6Class(
  classname = 'ConfigList',
  inherit = MappingBase,
  public = list(
    set = function(...) {
      args <- list(...)
      stopifnot(is_named_list(x = args))
      for (i in seq_along(along.with = args)) {
        key <- names(x = args)[i]
        value <- args[[i]]
        if (inherits(x = value, what = 'ScalarMap')) {
          self$add_map(key = key, value = value)
        } else {
          self$set_from_list(key = key, key = as.list(x = value))
        }
        return(invisible(x = self))
      }
      # stopifnot(is.character(x = key), length(x = key) == 1L)
      # if (is.null(x = value)) {
      #   super$set(key = key, value = value)
      #   return(invisible(x = self))
      # }
      # if (inherits(x = value, what = 'ScalarMap')) {
      #   self$add_map(key = key, value = value)
      #   return(invisible(x = self))
      # }
      # self$set_from_list(key = key, value = as.list(x = value))
    },
    add_map = function(key, value) {
      stopifnot(
        is.character(x = key),
        length(x = key) == 1L,
        inherits(x = value, what = 'ScalarMap')
      )
      private$.data[[key]] <- value
      return(invisible(x = self))
    },
    set_from_list = function(key, value) {
      stopifnot(is_named_list(x = value))
      cfg <- ScalarMap$new()
      for (i in names(x = value)) {
        cfg$set(key = i, value = value[[i]])
      }
      self$add_map(key = key, value = cfg)
      return(invisible(x = self))
    }
  )
)

setClass(Class = 'ConfigListS4', contains = 'list')

#' @method [ ConfigListS4
#' @export
#'
'[.ConfigListS4' <- function(x, i, ...) {
  .NotYetImplemented()
}

#' @method [[ ConfigListS4
#' @export
#'
'[[.ConfigListS4' <- function(x, i, ...) {
  .NotYetImplemented()
}

setAs(
  from = 'list',
  to = 'ConfigListS4',
  def = function(from) {
    stopifnot(!is.object(x = from))
  }
)

setMethod(
  f = '[[<-',
  signature = c(
    x = 'ConfigListS4',
    i = 'character',
    j = 'missing',
    value = 'ANY'
  ),
  definition = function(x, i, ..., value) {
    .NotYetImplemented()
  }
)

setMethod(
  f = '[[<-',
  signature = c(
    x = 'ConfigListS4',
    i = 'character',
    j = 'missing',
    value = 'data.frame'
  ),
  definition = function(x, i, ..., value) {
    .NotYetImplemented()
  }
)

setMethod(
  f = '[[<-',
  signature = c(
    x = 'ConfigListS4',
    i = 'character',
    j = 'missing',
    value = 'list'
  ),
  definition = function(x, i, ..., value) {
    .NotYetImplemented()
  }
)

setValidity(
  Class = 'ConfigListS4',
  method = function(object) {
    valid <- NULL
    return(valid %||% TRUE)
  }
)
