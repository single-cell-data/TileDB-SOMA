#' R6 Base Mapping Type
#'
#' Virtual base mapping type for R6 objects; defines internal data structure
#' (\code{private$.data}) as a named list along with behavior methods for
#' getting (\code{self$get()}) and setting (\code{self$set()}) items in the map
#'
#' @keywords internal
#' @export

MappingBase <- R6::R6Class(
  classname = 'MappingBase',
  lock_class = TRUE,
  public = list(
    #' @param ... Ignored
    #'
    #' @return NOPENOPENOPE[results=rd]{tiledbsoma:::rd_return_virtual()}
    #'
    initialize = function(...) {
      calls <- vapply_char(
        X = lapply(X = sys.calls(), FUN = as.character),
        FUN = '[[',
        1L
      )
      if ('MappingBase$new' %in% calls) {
        stop(
          "'MappingBase' is a virtual class and cannot be instantiated directly",
          call. = FALSE
        )
      }
      private$.data <- list()
    },
    #' @return The keys of the map
    #'
    keys = function() {
      return(names(private$.data))
    },
    #' @return A `list` containing the map values
    #'
    values = function() {
      return(unname(self$items()))
    },
    #' @return Return the items of the map as a list
    #'
    items = function() {
      return(private$.data)
    },
    #' @param key Key to fetch
    #' @templateVar key key
    #' @templateVar default NULL
    #' @template param-default
    #'
    #' @return The value of \code{key} in the map, or \code{default} if
    #' \code{key} is not found
    #'
    get = function(key, default = quote(expr = )) {
      stopifnot(
        "'key' must be a single character value" = is_scalar_character(key)
      )
      value <- private$.data[[key]]
      if (is.null(value)) {
        if (missing(default) || identical(x = default, y = quote(expr = ))) {
          private$.key_error(key)
        }
        return(default)
      }
      return(value)
    },
    #' @param key Key to set
    #' @templateVar key key
    #' @template param-value
    #'
    #' @return \[chainable\] Invisibly returns \code{self} with
    #' \code{value} added as \code{key}
    #'
    set = function(key, value) {
      stopifnot(
        "'key' must be a single character value" = is_scalar_character(key)
      )
      private$.data[[key]] <- value
      private$.data <- Filter(f = length, x = private$.data)
      if (!self$length()) {
        private$.data <- unname(obj = private$.data)
      }
      return(invisible(self))
    },
    #' @param ... Named arguments to add to \code{self}
    #'
    #' @return \[chainable\] Invisibly returns \code{self} with the values
    #' of \code{...} added to the map
    #'
    setv = function(...) {
      args <- as.list(c(...))
      stopifnot("all arguments must be named" = is_named_list(args))
      for (i in seq_along(args)) {
        self$set(key = names(args)[i], value = args[[i]])
      }
      return(invisible(self))
    },
    #' @param key Key to remove
    #'
    #' @return \[chainable\] Invisibly returns \code{self} with \code{key}
    #' removed from the map
    #'
    remove = function(key) {
      self$set(key = key, value = NULL)
      return(invisible(self))
    },
    #' @param map A mapping type to update the current map with
    #'
    #' @return \[chainable\] Invisibly returns \code{self} with the value
    #' of \code{map}
    #'
    update = function(map) {
      stopifnot(
        "'map' must be a mapping type" = inherits(x = map, what = 'MappingBase')
      )
      self$setv(map$items())
      return(invisible(self))
    },
    #' @return The number of items in the map
    #'
    length = function() {
      return(length(private$.data))
    },
    #' @return The map as a list
    #'
    to_list = function() {
      return(self$items())
    },
    #' @return \[chainable\] Prints information about the map to \code{stdout}
    #' and invisibly returns \code{self}
    #'
    print = function() {
      cat("<", class(self)[1L], ">\n", sep = '')
      if (length(self)) {
        cat(
          '  ',
          paste(self$keys(), self$values(), sep = ': ', collapse = '\n  '),
          '\n',
          sep = ''
        )
      }
      return(invisible(x = self))
    }
  ),
  private = list(
    .data = list(),
    .key_error = function(key) {
      stop("No key named ", sQuote(key), " found", call. = FALSE)
    }
  )
)

#' @param x A mapping object
#' @param i A key to fetch or set; see \code{$get()} or \code{$set()}
#' methods below
#' @template param-dots-ignored
#' @templateVar key i
#' @templateVar default NULL
#' @template param-default
#'
#' @return \code{[[}: The value of \code{i} in the map, or \code{default} if
#' \code{i} is not found
#'
#' @rdname MappingBase
#'
#' @method [[ MappingBase
#' @export
#'
'[[.MappingBase' <- function(x, i, ..., default = NULL) {
  return(x$get(key = i, default = default))
}

#' @templateVar key i
#' @template param-value
#'
#' @return \code{[[<-}: \code{x} with \code{value} added as \code{i}
#'
#' @rdname MappingBase
#'
#' @method [[<- MappingBase
#' @export
#'
'[[<-.MappingBase' <- function(x, i, ..., value) {
  stopifnot("'i' must be a single character value" = is_scalar_character(i))
  x$set(key = i, value = value)
  return(x)
}

#' @return \code{as.list}: The map as a list
#'
#' @rdname MappingBase
#'
#' @method as.list MappingBase
#' @export
#'
as.list.MappingBase <- function(x, ...) {
  return(x$to_list())
}

#' @return \code{length}: The number of items in the map
#'
#' @rdname MappingBase
#'
#' @method length MappingBase
#' @export
#'
length.MappingBase <- function(x) {
  return(x$length())
}

#' @return \code{names}: The keys of the map
#'
#' @rdname MappingBase
#'
#' @method names MappingBase
#' @export
#'
names.MappingBase <- function(x) {
  return(x$keys())
}
