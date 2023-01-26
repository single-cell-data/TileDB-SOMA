#' @include SOMAContextBase.R
#' @importClassesFrom tiledb tiledb_ctx
#'
NULL

#' @exportClass SOMATileDBContext
#'
setClass(
  Class = 'SOMATileDBContext',
  contains = c('SOMAContextBase', 'tiledb_ctx')
)

#' @method [[ SOMATileDBContext
#' @export
#'
'[[.SOMATileDBContext' <- function(x, i, ...) {
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

#' @method names SOMATileDBContext
#' @export
#'
names.SOMATileDBContext <- function(x) {
  nn <- NextMethod()
  return(c(nn, context_names(x = x)))
}

setMethod(
  f = '[[<-',
  signature = c(
    x = 'SOMATileDBContext',
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
  signature = 'SOMATileDBContext',
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
    if (length(x = soma_opts)) {
      .Object[names(x = soma_opts)] <- soma_opts
    }
    # Validate and return
    validObject(object = .Object)
    return(.Object)
  }
)

context_names <- function(x) {
  stopifnot(inherits(x = x, what = 'tiledb_ctx'))
  return(names(x = as.vector(x = tiledb::config(object = x))))
}
