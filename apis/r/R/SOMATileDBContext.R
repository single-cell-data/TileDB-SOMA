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
        !length(x = config) || is_named(x = config, allow_empty = FALSE)
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
      private$.tiledb_ctx <- tiledb::tiledb_ctx(config = cfg, cached = cached)
    },
    #' @return The keys of the map
    #'
    keys = function() {
      return(c(super$keys(), private$.tiledb_ctx_names()))
    },
    #' @return Return the items of the map as a list
    #'
    items = function() {
      return(c(super$items(), as.list(x = tiledb::config(object = private$.tiledb_ctx))))
    },
    #' @return The number of items in the map
    #'
    length = function() {
      return(super$length() + length(x = private$.tiledb_ctx_names()))
    },
    #' @param key Key to fetch
    #' @templateVar key key
    #' @templateVar default NULL
    #' @template param-default
    #'
    #' @return The value of \code{key} in the map, or \code{default} if
    #' \code{key} is not found
    #'
    #' @importFrom rlang missing_arg is_missing
    get = function(key, default = rlang::missing_arg()) {
      key <- key[1L]
      key <- tryCatch(
        expr = match.arg(arg = key, choices = self$keys()),
        error = \(...) NULL
      )
      if (is.null(x = key)) {
        if (rlang::is_missing(x = default)) {
          private$.key_error(key = key)
        }
      }
      if (key %in% private$.tiledb_ctx_names()) {
        val <- tiledb::config(object = private$.tiledb_ctx)[key]
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
        "'key' must be a single character value" = is_scalar_character(key),
        "'value' must be a single atomic value" = is.atomic(x = value) && length(value) == 1L,
      )
      if (key %in% private$.tiledb_ctx_names()) {
        cfg <- tiledb::config(object = private$.tiledb_ctx)
        cfg[key] <- as.character(x = value)
        private$.tiledb_ctx <- tiledb::tiledb_ctx(config = cfg)
      } else {
        super$set(key = key, value = value)
      }
      return(invisible(x = self))
    },
    #' @return A \code{\link[tiledb:tiledb_ctx]{tiledb_ctx}} object, dynamically
    #' constructed. Most useful for the constructor of this class.
    #'
    to_tiledb_context = function() {
      items <- sapply(X = super$items(), FUN = as.character, USE.NAMES = TRUE)
      cfg <- tiledb::config(object = private$.tiledb_ctx)
      for (opt in names(x = items)) {
        cfg[opt] <- items[opt]
      }
      tiledb_ctx <- tiledb::tiledb_ctx(config = cfg)
      return(tiledb_ctx)
    },
    #' @return A \code{\link[tiledb:tiledb_ctx]{tiledb_ctx}} object, which is
    #' a stored (and long-lived) result from \code{to_tiledb_context}.
    get_tiledb_context = function() {
      return(private$.tiledb_ctx)
    }
  ),
  private = list(
    .tiledb_ctx = NULL,
    .tiledb_ctx_names = function() {
      if (!inherits(x = private$.tiledb_ctx, what = 'tiledb_ctx')) {
        return(NULL)
      }
      return(tryCatch(
        expr = names(x = as.vector(x = tiledb::config(object = private$.tiledb_ctx))),
        error = \(...) NULL
      ))
    }
  )
)
