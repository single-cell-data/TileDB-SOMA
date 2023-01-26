#' @importFrom methods callNextMethod initialize setClass setMethod
#' setValidity slot slot<-
#'
NULL

#' @exportClass SOMAContextBase
#'
setClass(Class = 'SOMAContextBase', contains = c('VIRTUAL', 'list'))

#' @method [ SOMAContextBase
#' @export
#'
'[.SOMAContextBase' <- function(x, i, ...) {
  .NotYetImplemented()
}

#' @method [[ SOMAContextBase
#' @export
#'
"[[.SOMAContextBase" <- function(x, i, ...) {
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
  signature = c(x = 'SOMAContextBase'),
  definition = function(x, i, ..., value) {
    .NotYetImplemented()
  }
)

setMethod(
  f = '[[<-',
  signature = c(x = 'SOMAContextBase', i = 'ANY', j = 'missing', value = 'ANY'),
  definition = function(x, i, ..., value) {
    .NotYetImplemented()
  }
)

setMethod(
  f = '[[<-',
  signature = c(
    x = 'SOMAContextBase',
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
  Class = 'SOMAContextBase',
  method = function(object) {
    valid <- NULL
    # Ensure all entries are named
    if (is.null(x = names(x = object)) || !all(nzchar(x = names(x = object)))) {
      valid <- c(valid, "All entries must be named")
    }
    soma_contexts <- .SOMA_CONTEXTS()
    # Check individual entries
    for (i in names(x = object)) {
      val <- slot(object = object, name = '.Data')[[i]]
      if (length(x = val) != 1L) {
        valid <- c(
          valid,
          paste0(
            "All entries must be of length 1 (offending: ",
            sQuote(x = i),
            ")"
          )
        )
      }
      if (!is.atomic(x = val)) {
        valid <- c(
          valid,
          paste0(
            "All entries must be atomic (",
            sQuote(x = i),
            "is of type ",
            class(x = val)[1L],
            ")"
          )
        )
      }
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
