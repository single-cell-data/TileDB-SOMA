#' @include SOMAContextBase.R
#' @importClassesFrom tiledb tiledb_ctx
#'
NULL

SOMATileDBContext <- R6::R6Class(
  classname = 'SOMATileDBContext',
  inherit = SOMAContextBase,
  public = list(
    initialize = function(config = NULL, cached = TRUE) {
      super$initialize()
      config <- config %||% character()
      config['sm.mem.reader.sparse_global_order.ratio_array_data'] <- '0.3'
      stopifnot(is.character(x = config), !is.null(x = names(x = config)))
      # Identify options that are SOMA-specific
      soma_opts <- which(x = names(x = config) %in% names(x = .SOMA_CONTEXTS()))
      if (length(x = soma_opts)) {
        soma_config <- config[soma_opts]
        config <- config[-soma_opts]
      } else {
        soma_config <- vector()
      }
      # Add the TileDB context
      if (!length(x = config)) {
        config <- NA_character_
      }
      cfg <- tiledb::tiledb_config(config = config)
      private$.ctx <- tiledb::tiledb_ctx(config = cfg, cached = cached)
      # Add the SOMA options
      if (length(x = soma_config)) {
        self$setv(soma_config)
      }
    },
    keys = function() {
      return(super$keys(), private$.ctx_names())
    },
    get = function(key, default = NULL) {
      key <- match.arg(arg = key, choices = self$keys())
      if (key %in% private$.ctx_names()) {
        val <- tiledb::config(object = private$.ctx)[key]
        names(x = val) <- key
        return(val)
      }
      return(super$get(key = key, default = default))
    },
    set = function(key, value) {
      stopifnot(
        is.character(x = key),
        is.atomic(x = value),
        length(x = key) == length(x = value),
        length(x = key) == 1L
      )
      if (key %in% private$.ctx_names()) {
        cfg <- tiledb::config(object = private$.ctx)
        cfg[key] <- as.character(x = value)
        private$.ctx <- tiledb::tiledb_ctx(config = cfg)
      } else {
        super$set(key = key, value = value)
      }
      return(invisible(x = self))
    }
  ),
  private = list(
    .ctx = NULL,
    .ctx_names = function() {
      if (!inherits(x = private$.ctx, what = 'tiledb_ctx')) {
        return(NULL)
      }
      return(tryCatch(
        expr = names(x = as.vector(x = tiledb::config(object = private$.ctx))),
        error = \(...) NULL
      ))
    }
  )
)

setClass(
  Class = 'SOMATileDBContextS4',
  contains = c('SOMAContextBaseS4', 'tiledb_ctx')
)

#' @method [[ SOMATileDBContextS4
#' @export
#'
'[[.SOMATileDBContextS4' <- function(x, i, ...) {
  i <- i[1L]
  i <- tryCatch(
    expr = match.arg(arg = i, choices = names(x = x)),
    error = \(...) i
  )
  if (i %in% context_names(x = x)) {
    val <- tiledb::config(object = x)[i]
    names(x = val) <- i
    return(val)
  }
  return(NextMethod())
}

#' @method names SOMATileDBContextS4
#' @export
#'
names.SOMATileDBContextS4 <- function(x) {
  return(c(NextMethod(), context_names(x = x)))
}

setMethod(
  f = '[[<-',
  signature = c(
    x = 'SOMATileDBContextS4',
    i = 'character',
    j = 'missing',
    value = 'ANY'
  ),
  definition = function(x, i, ..., value) {
    stopifnot(length(x = i) == length(x = value), is.atomic(x = value))
    for (idx in seq_along(along.with = i)) {
      ii <- i[idx]
      vi <- value[idx]
      if (ii %in% context_names(x = x)) {
        cfg <- tiledb::config(object = x)
        cfg[ii] <- vi
        slot(object = x, name = 'ptr') <- slot(
          object = tiledb::tiledb_ctx(config = cfg),
          name = 'ptr'
        )
      } else {
        x <- callNextMethod(x = x, i = ii, value = vi)
      }
    }
    validObject(object = x)
    return(x)
  }
)

setMethod(
  f = 'initialize',
  signature = 'SOMATileDBContextS4',
  definition = function(.Object, config = NULL, cached = TRUE, ...) {
    .Object <- callNextMethod(.Object, ...)
    # tiledb::tiledb_config requires characters only?
    config <- config %||% character()
    config['sm.mem.reader.sparse_global_order.ratio_array_data'] <- '0.3'
    stopifnot(is.character(x = config) && !is.null(x = names(x = config)))
    # Identify options that are SOMA-specific
    soma_opts <- which(x = names(x = config) %in% names(x = .SOMA_CONTEXTS()))
    if (length(x = soma_opts)) {
      soma_config <- config[soma_opts]
      config <- config[-soma_opts]
    } else {
      soma_config <- c()
    }
    # Add the TileDB context
    if (!length(x = config)) {
      config <- NA_character_
    }
    cfg <- tiledb::tiledb_config(config = config)
    ctx <- tiledb::tiledb_ctx(config = cfg, cached = cached)
    slot(object = .Object, name = 'ptr') <- slot(object = ctx, name = 'ptr')
    # Add the SOMA options
    slot(object = .Object, name = '.Data') <- list()
    if (length(x = soma_config)) {
      .Object[names(x = soma_config)] <- soma_config
    }
    # Validate and return
    validObject(object = .Object)
    return(.Object)
  }
)

context_names <- function(x) {
  stopifnot(inherits(x = x, what = 'tiledb_ctx'))
  return(tryCatch(
    expr = names(x = as.vector(x = tiledb::config(object = x))),
    error = \(...) NULL
  ))
}
