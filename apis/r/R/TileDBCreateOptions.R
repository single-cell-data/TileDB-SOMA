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

# Non-filter-related schema parameters
DEFAULT_TILE_ORDER  <- function() { "ROW_MAJOR" }
DEFAULT_CELL_ORDER  <- function() { "ROW_MAJOR" }
DEFAULT_TILE_EXTENT <- function() { 2048 }
DEFAULT_CAPACITY    <- function() { 100000 }

# Filter-related schema parameters
DEFAULT_DATAFRAME_DIM_ZSTD_LEVEL       <- function() { 3 }
DEFAULT_SPARSE_ND_ARRAY_DIM_ZSTD_LEVEL <- function() { 3 }
DEFAULT_OFFSETS_FILTERS  <- function() { list("DOUBLE_DELTA", "BIT_WIDTH_REDUCTION", "ZSTD")}
DEFAULT_VALIDITY_FILTERS <- function() { list() }

# Used for chunked data ingestion
DEFAULT_WRITE_X_CHUNKED  <- function() { TRUE }
DEFAULT_GOAL_CHUNK_NNZ   <- function() { 200000000 }

#' TileDBCreateOptions
#'
#' Provides strongly-typed access and default values for `platform_config` options stored under the
#' "tiledb"->"create" mapping keys.
#'
#' Intended for internal use only.
#'
#' @noMd
#'
TileDBCreateOptions <- R6::R6Class(
  classname = 'TileDBCreateOptions',
  inherit = MappingBase,
  public = list(

    # Initializes from a PlatformConfig object.
    initialize = function(platform_config) {
      if (!is.null(platform_config)) {
        stopifnot(
          "'platform_config' must be null, or a PlatformConfig type" = inherits(x = platform_config, what = 'PlatformConfig')
        )

        if ("tiledb" %in% platform_config$platforms() && "create" %in% platform_config$params()) {
          map <- platform_config$get('tiledb', 'create')
          for (key in map$keys()) {
            super$set(key, map$get(key))
          }
        }
      }
    },


    #' @return int
    dataframe_dim_zstd_level = function() {
      super$get("dataframe_dim_zstd_level", default=DEFAULT_DATAFRAME_DIM_ZSTD_LEVEL())
    },

    #' @return int
    sparse_nd_array_dim_zstd_level = function() {
        super$get("sparse_nd_array_dim_zstd_level", default=DEFAULT_SPARSE_ND_ARRAY_DIM_ZSTD_LEVEL())
    },

    #' Returns the cell and tile orders that should be used.
    #' If neither ``cell_order`` nor ``tile_order`` is present, only in this
    #' case will we use the default values provided.
    #
    #' @return a list keyed by `"cell_order"` and `"tile_order"`, where either
    # value or both may be NULL.
    cell_tile_orders = function() {
        if ("cell_order" %in% super$keys() || "tile_order" %in% super$keys()) {
            c(cell_order = super$get("cell_order", default=NULL), tile_order = super$get("tile_order", default=NULL))
        } else {
          c(cell_order = DEFAULT_CELL_ORDER(), tile_order = DEFAULT_TILE_ORDER())
        }
    },

    # This allows a user to specify TileDB extents on a per-dimension basis.
    #
    # Example:
    # * cfg <- PlatformConfig$new()
    # * cfg$set('tiledb', 'create', 'dims', list(soma_dim_0 = list(tile = 999)))
    # * tdco <- TileDBCreateOptions$new(cfg)
    #
    # Then tdco$dim_tile("soma_dim_0") will be 999
    #
    #' @return int
    dim_tile = function(dim_name, default=DEFAULT_TILE_EXTENT()) {
      stopifnot(!is.null(dim_name))
      o <- self$.dim(dim_name)
      if (is.null(o)) {
        return(default)
      }
      o <- o[['tile']]
      if (is.null(o)) {
        return(default)
      }
      return(o)
    },

    #' @return int
    capacity = function() {
      super$get("capacity", default=DEFAULT_CAPACITY())
    },

    # Example input:
    #
    #   cfg <- PlatformConfig$new()
    #   cfg$set('tiledb', 'create', 'attrs', list(
    #     soma_dim_0 = list(tile = 100, filters = list("RLE")),
    #     soma_dim_1 = list(tile = 200, filters = list("RLE", list(name="ZSTD", COMPRESSION_LEVEL=9)))
    #   ))
    #   tdco <- TileDBCreateOptions$new(cfg)
    #
    #' @return list of tiledb.Filter
    dim_filters = function(dim_name, default=list()) {
      stopifnot(!is.null(dim_name))
      stopifnot(!is.null(default))
      o <- self$.dim(dim_name)
      if (is.null(o)) {
        return(self$.build_filters(default))
      }
      o <- o[['filters']]
      if (is.null(o)) {
        return(self$.build_filters(default))
      }
      return(self$.build_filters(o))
    },

    # Example input:
    #
    #   cfg <- PlatformConfig$new()
    #   cfg$set('tiledb', 'create', 'attrs', list(
    #     soma_data_a = list(filters = list("RLE")),
    #     soma_data_b = list(filters = list("RLE", list(name="ZSTD", COMPRESSION_LEVEL=9)))
    #   ))
    #   tdco <- TileDBCreateOptions$new(cfg)
    #
    #' @return list of tiledb.Filter
    attr_filters = function(attr_name, default=list()) {
      stopifnot(!is.null(attr_name))
      stopifnot(!is.null(default))
      o <- self$.attr(attr_name)
      if (is.null(o)) {
        return(self$.build_filters(default))
      }
      o <- o[['filters']]
      if (is.null(o)) {
        return(self$.build_filters(default))
      }
      return(self$.build_filters(o))
    },

    #' @return list of tiledb.Filter
    offsets_filters = function(default=DEFAULT_OFFSETS_FILTERS()) {
      self$.build_filters(super$get("offsets_filters", default))
    },

    #' @return list of tiledb.Filter
    validity_filters = function(default=DEFAULT_VALIDITY_FILTERS()) {
      self$.build_filters(super$get("validity_filters", default))
    },

    #' @return bool
    write_X_chunked = function() {
        super$get("write_X_chunked", default=DEFAULT_WRITE_X_CHUNKED())
    },

    #' @return int
    goal_chunk_nnz = function() {
        super$get("goal_chunk_nnz", default=DEFAULT_GOAL_CHUNK_NNZ())
    },

    # This is an accessor for nested things like
    #
    #   cfg <- PlatformConfig$new()
    #   cfg$set('tiledb', 'create', 'dims', list(soma_dim_0 = list(tile = 6)))
    #   tdco <- TileDBCreateOptions$new(cfg)
    #
    # -- given "soma_dim_0", this pulls out the `list(tile = 6)` part.
    #
    #' return Named list of character
    .dim = function(dim_name) {
        o <- super$get("dims", NULL)
        if (is.null(o)) {
          return(NULL)
        }
        o[[dim_name]]
    },

    # This is an accessor for nested things like
    #
    #   cfg <- PlatformConfig$new()
    #   cfg$set('tiledb', 'create', 'attrs', list(myattr = list(foo = 'bar')))
    #   tdco <- TileDBCreateOptions$new(cfg)
    #
    # -- this pulls out the `attrs` -> `myattr` part.
    #
    #' return Named list of character
    .attr = function(attr_name) {
        o <- super$get("attrs", NULL)
        if (is.null(o)) {
          return(NULL)
        }
        o[[attr_name]]
    },

    # The `items` argument is a list of arguments acceptable to `.build_filter`.
    # Example argument:
    #   list(
    #     "RLE",
    #     list(name = "ZSTD", COMPRESSION_LEVEL = 9)
    #   )
    # Returns a list of tiledb filters.
    .build_filters = function(items) {
      lapply(items, self$.build_filter)
    },

    # The `item` argument can be a string, like "RLE".
    # Or, a named list with filter name and remaining arguments, like
    # list(name="ZSTD", COMPRESSION_LEVEL=-1).
    #
    # Returns a tiledb filter.
    .build_filter = function(item) {
      # See also:
      # https://tiledb-inc.github.io/TileDB-R/reference/tiledb_filter.html
      if (is.character(item) && length(item) == 1) {
        return(tiledb::tiledb_filter(item[[1]]))
      }

      stopifnot(!is.na(item['name']))
      filter <- tiledb::tiledb_filter(item[['name']])
      for (key in names(item)) {
        if (key != "name") {
          tiledb::tiledb_filter_set_option(filter, key, item[[key]])
        }
      }
      filter
    }
  )
)
