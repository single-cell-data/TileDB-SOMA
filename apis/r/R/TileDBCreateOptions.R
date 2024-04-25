# Example config (written JSON-style for compact reference -- see unit-test files for precise
# syntax):
#
# platform_config = {
#     "tiledb": {
#         "create": {
#             "cell_order": "ROW_MAJOR",
#             "tile_order": "COL_MAJOR",
#             "capacity": 8888,
#             "dataframe_dim_zstd_level": 1,
#             "sparse_nd_array_dim_zstd_level": 2,
#             "offsets_filters": ["RLE", "NONE"],
#             "dims": {
#                 "soma_dim_0": {"tile": 6},
#                 # Empty filters for soma_dim_1 overrides the default
#                 # dimension zstd level defined below.
#                 "soma_dim_1": {"filters": []},
#             },
#             "attrs": {"soma_data": {"filters": ["RLE"]}},
#         }
#     }
# }

.CREATE_DEFAULTS <- list(
  # Non-filter-related schema parameters
  tile_order = 'ROW_MAJOR',
  cell_order = 'ROW_MAJOR',
  # tile_extent = 2048,
  capacity = 100000,
  allows_duplicates = FALSE,
  # Filter-related schema parameters
  dataframe_dim_zstd_level = 3,
  sparse_nd_array_dim_zstd_level = 3,
  offsets_filters = list('DOUBLE_DELTA', 'BIT_WIDTH_REDUCTION', 'ZSTD'),
  validity_filters = list(),
  # Used for chunked data ingestion
  write_X_chunked = TRUE,
  goal_chunk_nnz = 200000000
)

#' TileDBCreateOptions
#'
#' @description Provides strongly-typed access and default values for
#' \code{platform_config} options stored under the \dQuote{tiledb}
#' \eqn{\rightarrow} \dQuote{create} mapping keys.
#'
#' Intended for internal use only.
#'
#' @keywords internal
#'
#' @export
#'
#' @noMd
#'
TileDBCreateOptions <- R6::R6Class(
  classname = 'TileDBCreateOptions',
  inherit = MappingBase,
  public = list(

    #' @description Create a \code{TileDBCreateOptions} object
    #'
    #' @template param-platform-config
    #'
    initialize = function(platform_config = NULL) {
      if (!is.null(platform_config)) {
        stopifnot("'platform_config' must be a PlatformConfig" = inherits(
          x = platform_config,
          what = 'PlatformConfig'
        ))
        if ("tiledb" %in% platform_config$platforms() && "create" %in% platform_config$params()) {
          map <- platform_config$get('tiledb', 'create')
          for (key in map$keys()) {
            super$set(key, map$get(key))
          }
        }
      }

      for (key in setdiff(names(.CREATE_DEFAULTS), self$keys())) {
        self$set(key = key, value = .CREATE_DEFAULTS[[key]])
      }
    },

    #' @description Returns the cell and tile orders that should be used.
    #' If neither \code{cell_order} nor \code{tile_order} is present, only in this
    #' case will we use the default values provided.
    #
    #' @return A two-length character vector with names of
    #' \dQuote{\code{cell_order}} and \dQuote{\code{tile_order}}
    #'
    cell_tile_orders = function() c(
      cell_order = self$get('cell_order'),
      tile_order = self$get('tile_order')
    ),

    #' @param dim_name Name of dimension to get tiling for
    #' @param default Default tiling if \code{dim_name} is not set
    #'
    #' @return int
    #'
    #' @examples
    #' cfg <- PlatformConfig$new()
    #' cfg$set(
    #'   platform = 'tiledb',
    #'   param = 'create',
    #'   key = 'dims',
    #'   value = list(soma_dim_0 = list(tile = 999))
    #' )
    #' (tdco <- TileDBCreateOptions$new(cfg))
    #' tdco$dim_tile("soma_dim_0")
    #' tdco$dim_tile("soma_dim_1")
    #'
    dim_tile = function(dim_name, default = 2048) {
      stopifnot(
        "'dim_name' must be a single character value" = is.character(dim_name) &&
          length(dim_name) == 1L &&
          nzchar(dim_name),
        "'default' must be a single integerish value" = rlang::is_scalar_integerish(
          default,
          finite = TRUE
        )
      )
      return(private$.dim(dim_name)[['tile']] %||% default)
    },

    #' @return int
    #'
    capacity = function() self$get('capacity'),

    #' @return bool
    #'
    allows_duplicates = function() self$get('allows_duplicates'),

    #' @return int
    #'
    dataframe_dim_zstd_level = function() self$get('dataframe_dim_zstd_level'),

    #' @return int
    #'
    sparse_nd_array_dim_zstd_level = function() self$get('sparse_nd_array_dim_zstd_level'),

    #' @param default Default offset filters to use if not currently set
    #'
    #' @return A list of
    #' \code{\link[tiledb:tiledb_filter-class]{tiledb_filter}} objects
    #'
    offsets_filters = function(default = list()) {
      stopifnot(
        "'default' must be an unnamed list" = is.list(default) && !is_named(default)
      )
      return(private$.build_filters(self$get('offsets_filters', default = default)))
    },

    #' @param default Default validity filters to use if not currently set
    #'
    #' @return A list of
    #' \code{\link[tiledb:tiledb_filter-class]{tiledb_filter}} objects
    #'
    validity_filters = function(default = list()) {
      stopifnot(
        "'default' must be an unnamed list" = is.list(default) && !is_named(default)
      )
      return(private$.build_filters(self$get('validity_filters', default = default)))
    },

    #' @param dim_name Name of dimension to get filters for
    #' @param default Default filters to use for if not currently set
    #'
    #' @return A list of
    #' \code{\link[tiledb:tiledb_filter-class]{tiledb_filter}} objects
    #'
    #' @examples
    #' filters <- list(
    #'   soma_dim_0 = list(tile = 100, filters = list("RLE")),
    #'   soma_dim_1 = list(tile = 200, filters = list("RLE", list(name = "ZSTD", COMPRESSION_LEVEL = 9)))
    #' )
    #' cfg <- PlatformConfig$new()
    #' cfg$set(platform = 'tiledb', param = 'create', key = 'dims', value = filters)
    #' (tdco <- TileDBCreateOptions$new(cfg))
    #' tdco$dim_filters("soma_dim_0")
    #' tdco$dim_filters("non-existant")
    #'
    dim_filters = function(dim_name, default = list()) {
      stopifnot(
        "'dim_name' must be a single character value" = is.character(dim_name) &&
          length(dim_name) == 1L &&
          nzchar(dim_name),
        "'default' must be an unnamed list" = is.list(default) && !is_named(default)
      )
      filters <- private$.dim(dim_name)[['filters']] %||% default
      return(private$.build_filters(filters))
    },

    #' @param attr_name Name of attribute
    #' @param default Default filters to use if not currently set
    #'
    #' @return A list of
    #' \code{\link[tiledb:tiledb_filter-class]{tiledb_filter}} objects
    #'
    #' @examples
    #' filters <- list(
    #'   soma_data_a = list(filters = list("RLE")),
    #'   soma_data_b = list(filters = list("RLE", list(name = "ZSTD", COMPRESSION_LEVEL = 9)))
    #' )
    #' cfg <- PlatformConfig$new()
    #' cfg$set(platform = 'tiledb', param = 'create', key = 'attrs', value = filters)
    #' (tdco <- TileDBCreateOptions$new(cfg))
    #' tdco$attr_filters("soma_data_b")
    #' tdco$attr_filters("non-existant")
    #'
    attr_filters = function(attr_name, default = list()) {
      stopifnot(
        "'attr_name' must be a single character value" = is.character(attr_name) &&
          length(attr_name) == 1L &&
          nzchar(attr_name),
        "'default' must be an unnamed list" = is.list(default) && !is_named(default)
      )
      filters <- private$.attr(attr_name)[['filters']] %||% default
      return(private$.build_filters(filters))
    },

    #' @return bool
    #'
    write_X_chunked = function() self$get('write_X_chunked'),

    #' @return int
    #'
    goal_chunk_nnz = function() self$get('goal_chunk_nnz'),

    #' @description ...
    #'
    #' @param build_filters Build filters into
    #' \code{\link[tiledb:tiledb_filter-class]{tiledb_filter}} objects
    #'
    #' @return The create options as a list
    #'
    to_list = function(build_filters = TRUE) {
      stopifnot("'build_filters' must be TRUE or FALSE" = is_scalar_logical(build_filters))
      opts <- super$to_list()
      for (key in grep('_filters$', names(.CREATE_DEFAULTS), value = TRUE)) {
        if (is.null(opts[[key]])) {
          opts[[key]] <- .CREATE_DEFAULTS[[key]]
        }
      }
      if (isTRUE(build_filters)) {
        for (key in grep('_filters$', names(x = opts), value = TRUE)) {
          opts[[key]] <- private$.build_filters(opts[[key]])
        }
        for (key in c('dims', 'attrs')) {
          for (i in seq_along(opts[[key]])) {
            if ('filters' %in% names(opts[[key]][[i]])) {
              opts[[key]][[i]][['filters']] <- private$.build_filters(
                opts[[key]][[i]][['filters']]
              )
            }
          }
        }
      }
      return(opts)
    }
  ),

  private = list(
    # This is an accessor for nested things like
    #
    #   cfg <- PlatformConfig$new()
    #   cfg$set('tiledb', 'create', 'dims', list(soma_dim_0 = list(tile = 6)))
    #   tdco <- TileDBCreateOptions$new(cfg)
    #
    # -- given "soma_dim_0", this pulls out the `list(tile = 6)` part.
    #
    # @param dim_name Name of dimensions
    #
    # @return Named list of character
    #
    .dim = function(dim_name) self$get('dims', NULL)[[dim_name]],

    # This is an accessor for nested things like
    #
    #   cfg <- PlatformConfig$new()
    #   cfg$set('tiledb', 'create', 'attrs', list(myattr = list(foo = 'bar')))
    #   tdco <- TileDBCreateOptions$new(cfg)
    #
    # -- this pulls out the `attrs` -> `myattr` part.
    #
    # @param attr_name Name of attribute
    #
    # @return Named list of character
    #'
    .attr = function(attr_name) self$get("attrs", NULL)[[attr_name]],

    # The `items` argument is a list of arguments acceptable to `.build_filter`.
    # Example argument:
    #   list(
    #     "RLE",
    #     list(name = "ZSTD", COMPRESSION_LEVEL = 9)
    #   )
    # @param items A list of filters; see \code{$.build_filter()} for details
    # about entries in \code{items}
    #
    # @return A list of tiledb filters.
    #
    .build_filters = function(items) lapply(items, private$.build_filter),

    # The `item` argument can be a string, like "RLE".
    # Or, a named list with filter name and remaining arguments, like
    # list(name="ZSTD", COMPRESSION_LEVEL=-1).
    #
    # @param item THe name of a filter or a list with the name and arguments
    # for a filter (eg. \code{list(name = "ZSTD", COMPRESSION_LEVEL = -1)})
    #
    # @return A tiledb filter.
    #
    .build_filter = function(item) {
      # See also:
      # https://tiledb-inc.github.io/TileDB-R/reference/tiledb_filter.html
      if (rlang::is_scalar_character(item)) {
        # return(tiledb::tiledb_filter(item[[1]]))
        item <- list(name = item)
      }
      stopifnot(
        "'item' must be a named list" = is.list(item) &&
          is_named(item, allow_empty = FALSE),
        "'name' must be one of the names in 'item'" = 'name' %in% names(item)
      )
      filter <- tiledb::tiledb_filter(item[['name']])
      for (key in setdiff(x = names(item), y = 'name')) {
        tiledb::tiledb_filter_set_option(filter, option = key, value = item[[key]])
      }
      return(filter)
    }
  )
)
