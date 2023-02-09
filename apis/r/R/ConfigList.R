#' A Configuration List
#'
#' @export
#'
ConfigList <- R6::R6Class(
  classname = 'ConfigList',
  inherit = MappingBase,
  public = list(
    #' @param op ...
    #' @param key ...
    #' @param value ...
    #'
    #' @return \[chainable\] ...
    set = function(op, key, value) {
      stopifnot(length(x = op) == 1L, is.character(x = op))
      opmap <- self$get(key = op, default = ScalarMap$new())
      if (missing(x = key) && inherits(x = value, what = 'ScalarMap')) {
        opmap$update(map = value)
        super$set(key = op, value = opmap)
        return(invisible(x = self))
      }
      stopifnot(length(x = key) == 1L, is.character(x = key))
      opmap$set(key = key, value = value)
      super$set(key = op, value = opmap)
      return(invisible(x = self))
    },
    setv = function(...) {
      .NotYetImplemented()
    }
    # add_map = function(key, value) {
    #   stopifnot(
    #     is.character(x = key),
    #     length(x = key) == 1L,
    #     inherits(x = value, what = 'ScalarMap')
    #   )
    #   private$.data[[key]] <- value
    #   return(invisible(x = self))
    # },
    # set_from_list = function(key, value) {
    #   stopifnot(is_named_list(x = value))
    #   cfg <- ScalarMap$new()
    #   for (i in names(x = value)) {
    #     cfg$set(key = i, value = value[[i]])
    #   }
    #   self$add_map(key = key, value = cfg)
    #   return(invisible(x = self))
    # }
  )
)
