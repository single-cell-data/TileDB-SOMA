#' @importFrom methods as callNextMethod initialize setAs setClass
#' setMethod setValidity slot slot<-
#'
#' @importFrom methods as setAs
#'
NULL

setClass(
  Class = 'ScalarList',
  contains = 'list',
  slots = list(type. = 'character')
)

#' @method [ ScalarList
#' @export
#'
'[.ScalarList' <- function(x, i, ...) {
  .NotYetImplemented()
}

setAs(
  from = 'list',
  to = 'ScalarList',
  def = function(from) {
    stopifnot(!is.object(x = from))
  }
)

setMethod(
  f = '[[<-',
  signature = c(
    x = 'ScalarList',
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
    x = 'ScalarList',
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
    x = 'ScalarList',
    i = 'character',
    j = 'missing',
    value = 'list'
  ),
  definition = function(x, i, ..., value) {
    .NotYetImplemented()
  }
)

setMethod(
  f = 'initialize',
  signature = 'ScalarList',
  definition = function(.Object, ..., type. = 'any') {
    .Object <- callNextMethod(.Object)
    # Set object type
    type. <- type.[1L]
    type. <- match.arg(arg = type., choices = .SCALAR_TYPES())
    slot(object = .Object, name = 'type.') <- type.
    # Add object data
    vals <- list(...)
    if (length(x = vals)) {
      stopifnot(is_named_list(vals))
    }
    # Coerce the data into specified type if needed
    if (type. != 'any') {
      vals <- sapply(X = vals, FUN = as, class = type., simplify = FALSE)
    }
    slot(object = .Object, name = '.Data') <- vals
    validObject(object = .Object)
    return(.Object)
  }
)

setValidity(
  Class = 'ScalarList',
  method = function(object) {
    valid <- NULL
    # Ensure bare object
    if (is.data.frame(x = object) || is.data.frame(x = slot(object = object, name = '.Data'))) {
      return("Object must be a bare list")
    }
    # Check type
    type. <- slot(object = object, name = 'type.')
    if (length(x = type.) != 1L) {
      valid <- c(valid, "'type.' must be a one-length character")
    } else if (!type. %in% .SCALAR_TYPES()) {
      valid <- c(valid, "'type.' must be a scalar type")
      is_check <- is.atomic
    } else {
      is_check <- switch(
        EXPR = type.,
        any = is.atomic,
        \(x) is.vector(x = x, mode = type.)
      )
    }
    # Check .Data
    if (length(x = object) && !is_named_list(x = object)) {
      valid <- c(valid, "Object must be a named list")
    }
    for (i in seq_len(length.out = length(x = object))) {
      ii <- names(x = object)[i]
      vi <- object[[i]]
      if (!is_check(x = vi) || length(x = vi) != 1L) {
        valid <- c(
          valid,
          paste0(
            "All entries must be scalars ",
            switch(EXPR = type., any = 'atomics', paste('of type', type.)),
            " (offending entry: ",
            ifelse(
              test = isTRUE(x = nzchar(x = ii)),
              yes = sQuote(x = ii),
              no = paste("#", i)
            ),
            ")"
          )
        )
      }
    }
    return(valid %||% TRUE)
  }
)

.SCALAR_TYPES <- function() {
  return(c('any', 'numeric', 'integer', 'character', 'logical'))
}
