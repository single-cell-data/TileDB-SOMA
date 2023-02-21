#' SOMA TileDB Context
#'
#' Context map for TileDB-backed SOMA objects
#'
#' @export
#'
#' @noMd
#'
SOMATileDBContext <- R6::R6Class(
  classname = 'SOMATileDBContext',
  inherit = SOMAContextBase,
  public = list(
    #' @template param-config
    #' @param cached Force new creation
    #'
    #' @return An instantiated \code{SOMATileDBContext} object
    #'
    initialize = function(config = NULL, cached = TRUE) {
      # super$initialize()
      config <- config %||% character()
      # Identify options that are SOMA-specific
      soma_opts <- which(x = names(x = config) %in% names(x = .SOMA_CONTEXTS()))
      if (length(x = soma_opts)) {
        soma_config <- config[soma_opts]
        config <- config[-soma_opts]
      } else {
        soma_config <- NULL
      }
      super$initialize(config = soma_config)
      if (is.list(x = config)) {
        config <- unlist(x = config)
      }
      stopifnot(
        !length(x = config) || is.character(x = config),
        !length(x = config) || is_named2(x = config)
      )
      config['sm.mem.reader.sparse_global_order.ratio_array_data'] <- '0.3'
      # Add the TileDB context
      # TODO: Reenable this with newer version of tiledb-r
      # if (!length(x = config)) {
      #   config <- NA_character_
      # }
      # cfg <- tiledb::tiledb_config(config = config)
      # TODO: replace this when patched version of tiledb-r is released
      cfg <- tiledb::tiledb_config()
      for (opt in names(x = config)) {
        cfg[opt] <- config[opt]
      }
      private$.ctx <- tiledb::tiledb_ctx(config = cfg, cached = cached)
    },
    #' @return The keys of the map
    #'
    keys = function() {
      return(c(super$keys(), private$.ctx_names()))
    },
    #' @return The number of items in the map
    #'
    length = function() {
      return(super$length() + length(x = private$.ctx_names()))
    },
    #' @param key Key to fetch
    #' @templateVar key key
    #' @templateVar default NULL
    #' @template param-default
    #'
    #' @return The value of \code{key} in the map, or \code{default} if
    #' \code{key} is not found
    #'
    get = function(key, default = NULL) {
      key <- match.arg(arg = key, choices = self$keys())
      if (key %in% private$.ctx_names()) {
        val <- tiledb::config(object = private$.ctx)[key]
        names(x = val) <- key
        return(val)
      }
      return(super$get(key = key, default = default))
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
