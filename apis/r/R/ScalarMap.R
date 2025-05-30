#' A Mapping Type for Scalars
#'
#' An R6 mapping type that is limited to scalar atomic vector types only; can
#' optionally be limited further to a specific atomic vector type
#' (eg. \dQuote{\code{logical}}).
#'
#' @keywords internal
#'
#' @export
#'
#' @examples
#' (map <- ScalarMap$new())
#' map$set("a", 1L)
#' map
#'
#' map$get("a")
#' map$get("b", default = NULL)
#'
ScalarMap <- R6::R6Class(
  classname = "ScalarMap",
  inherit = MappingBase,
  public = list(
    #' @param type Limit the \code{ScalarMap} to a preset type; choose from:
    #' \itemize{
    #'  \item \dQuote{\code{any}}
    #'  \item \dQuote{\code{numeric}}
    #'  \item \dQuote{\code{integer}}
    #'  \item \dQuote{\code{character}}
    #'  \item \dQuote{\code{logical}}
    #' }
    #'
    #' @return An instantiated \code{ScalarMap} object with the
    #' type set to \code{type}
    #'
    initialize = function(type = "any") {
      private$.type <- match.arg(arg = type, choices = .SCALAR_TYPES())
    },

    #' @param key Key to set
    #' @templateVar key key
    #' @template param-value
    #'
    #' @return \[chainable\] Invisibly returns \code{self} with
    #' \code{value} added as \code{key}
    #'
    set = function(key, value) {
      if (!is.null(value) && self$type != "any") {
        if (!inherits(x = value, what = self$type)) {
          stop(
            "'value' must be a ",
            sQuote(x = self$type, q = FALSE),
            call. = FALSE
          )
        }
      }
      super$set(key = key, value = value)
      return(invisible(self))
    }
  ),
  active = list(
    #' @field type The type that this \code{ScalarMap} is limited to
    #'
    type = function(value) {
      if (!missing(value)) {
        stop("The map type cannot be changed", call. = FALSE)
      }
      return(private$.type)
    }
  ),
  private = list(
    .type = character(1L)
  )
)

.SCALAR_TYPES <- function() {
  return(c("any", "numeric", "integer", "character", "logical"))
}
