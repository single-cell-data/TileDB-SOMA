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
.default_tile_order        <- function() { "ROW_MAJOR" }
.default_cell_order        <- function() { "ROW_MAJOR" }
.default_tile_extent       <- function() { 2048        }
.default_capacity          <- function() { 100000      }
.default_allows_duplicates <- function() { FALSE      }

# Filter-related schema parameters
.default_dataframe_dim_zstd_level       <- function() { 3 }
.default_sparse_nd_array_dim_zstd_level <- function() { 3 }
.default_offsets_filters  <- function() { list("DOUBLE_DELTA", "BIT_WIDTH_REDUCTION", "ZSTD")}
.default_validity_filters <- function() { list() }

# Used for chunked data ingestion
.default_write_x_chunked  <- function() { TRUE }
.default_goal_chunk_nnz   <- function() { 200000000 }

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

    # Initializes from a PlatformConfig object.
    #' @description Create a \code{TileDBCreateOptions} object
    #'
    #' @template param-platform-config
    #'
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

    #' @description Returns the cell and tile orders that should be used.
    #' If neither \code{cell_order} nor \code{tile_order} is present, only in this
    #' case will we use the default values provided.
    #
    #' @return a list keyed by \dQuote{\code{cell_order}} and
    #' \dQuote{\code{tile_order}}, where either value or both may be \code{NULL}.

    cell_tile_orders = function() {
        if ("cell_order" %in% super$keys() || "tile_order" %in% super$keys()) {
            c(cell_order = super$get("cell_order", default=NULL), tile_order = super$get("tile_order", default=NULL))
        } else {
          c(cell_order = .default_cell_order(), tile_order = .default_tile_order())
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
    #' @param dim_name Name of dimension to get tiling for
    #' @param default Default tiling if \code{dim_name} is not set
    #'
    #' @return int
    #'
    dim_tile = function(dim_name, default=.default_tile_extent()) {
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
      super$get("capacity", default=.default_capacity())
    },

    #' @return bool
    allows_duplicates = function() {
      super$get("allows_duplicates", default=.default_allows_duplicates())
    },

    #' @return int
    dataframe_dim_zstd_level = function() {
      super$get("dataframe_dim_zstd_level", default=.default_dataframe_dim_zstd_level())
    },

    #' @return int
    sparse_nd_array_dim_zstd_level = function() {
        super$get("sparse_nd_array_dim_zstd_level", default=.default_sparse_nd_array_dim_zstd_level())
    },

    #' @param default Default offset filters to use if not currently set
    #'
    #' @return list of tiledb.Filter
    offsets_filters = function(default=.default_offsets_filters()) {
      self$.build_filters(super$get("offsets_filters", default))
    },

    #' @param default Default validity filters to use if not currently set
    #'
    #' @return list of tiledb.Filter
    validity_filters = function(default=.default_validity_filters()) {
      self$.build_filters(super$get("validity_filters", default))
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
    #' @param dim_name Name of dimension to get filters for
    #' @param default Default filters to use for if not currently set
    #'
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
    #' @param attr_name Name of attribute
    #' @param default Default filters to use if not currently set
    #'
    #' @return list of tiledb.Filter
    #'
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

    #' @return bool
    #'
    write_X_chunked = function() {
        super$get("write_X_chunked", default=.default_write_x_chunked())
    },

    #' @return int
    #'
    goal_chunk_nnz = function() {
        super$get("goal_chunk_nnz", default=.default_goal_chunk_nnz())
    },

    # This is an accessor for nested things like
    #
    #   cfg <- PlatformConfig$new()
    #   cfg$set('tiledb', 'create', 'dims', list(soma_dim_0 = list(tile = 6)))
    #   tdco <- TileDBCreateOptions$new(cfg)
    #
    # -- given "soma_dim_0", this pulls out the `list(tile = 6)` part.
    #
    #' @param dim_name Name of dimensions
    #'
    #' @return Named list of character
    #'
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
    #' @param attr_name Name of attribute
    #'
    #' @return Named list of character
    #'
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
    #' @param items A list of filters; see \code{$.build_filter()} for details
    #' about entries in \code{items}
    #'
    #' @return A list of tiledb filters.
    #'
    .build_filters = function(items) {
      lapply(items, self$.build_filter)
    },

    # The `item` argument can be a string, like "RLE".
    # Or, a named list with filter name and remaining arguments, like
    # list(name="ZSTD", COMPRESSION_LEVEL=-1).
    #
    #' @param item THe name of a filter or a list with the name and arguments
    #' for a filter (eg. \code{list(name = "ZSTD", COMPRESSION_LEVEL = -1)})
    #'
    #' @return A tiledb filter.
    #'
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
