#' SOMA TileDB Context
#'
#' @export
#'
SOMATileDBContext <- R6::R6Class(
  classname = 'SOMATileDBContext',
  inherit = SOMAContextBase,
  public = list(
    #' @param config ...
    #' @param cached ...
    #'
    #' @return ...
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
    #' @return ...
    keys = function() {
      return(super$keys(), private$.ctx_names())
    },
    #' @param key ...
    #' @param default ...
    #'
    #' @return ...
    get = function(key, default = NULL) {
      key <- match.arg(arg = key, choices = self$keys())
      if (key %in% private$.ctx_names()) {
        val <- tiledb::config(object = private$.ctx)[key]
        names(x = val) <- key
        return(val)
      }
      return(super$get(key = key, default = default))
    },
    #' @param key ...
    #' @param value ...
    #'
    #' @return ...
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
