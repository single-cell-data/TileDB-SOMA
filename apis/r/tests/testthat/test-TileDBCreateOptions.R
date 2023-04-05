# Verify that the TileDBCreateOptions takes a NULL or a PlatformConfig
test_that("TileDBCreateOptions construction", {
  expect_error(TileDBCreateOptions$new(c("foo": "bar")))

  tdco <- TileDBCreateOptions$new(NULL)

  cfg <- PlatformConfig$new()
  tdco <- TileDBCreateOptions$new(cfg)
})

# Verify access of parameters across a range of nestedness of keys in the PlatformConfig.
test_that("TileDBCreateOptions access from PlatformConfig", {
  cfg <- PlatformConfig$new()
  tdco <- TileDBCreateOptions$new(cfg)
  expect_equal(tdco$dataframe_dim_zstd_level(), DEFAULT_DATAFRAME_DIM_ZSTD_LEVEL())

  cfg <- PlatformConfig$new()
  cfg$set('not_tiledb', 'not_create', 'not_dataframe_dim_zstd_level', 999)
  tdco <- TileDBCreateOptions$new(cfg)
  expect_equal(tdco$dataframe_dim_zstd_level(), DEFAULT_DATAFRAME_DIM_ZSTD_LEVEL())

  cfg <- PlatformConfig$new()
  cfg$set('tiledb', 'not_create', 'not_dataframe_dim_zstd_level', 999)
  tdco <- TileDBCreateOptions$new(cfg)
  expect_equal(tdco$dataframe_dim_zstd_level(), DEFAULT_DATAFRAME_DIM_ZSTD_LEVEL())

  cfg <- PlatformConfig$new()
  cfg$set('not_tiledb', 'create', 'not_dataframe_dim_zstd_level', 999)
  tdco <- TileDBCreateOptions$new(cfg)
  expect_equal(tdco$dataframe_dim_zstd_level(), DEFAULT_DATAFRAME_DIM_ZSTD_LEVEL())

  cfg <- PlatformConfig$new()
  cfg$set('tiledb', 'create', 'not_dataframe_dim_zstd_level', 999)
  tdco <- TileDBCreateOptions$new(cfg)
  expect_equal(tdco$dataframe_dim_zstd_level(), DEFAULT_DATAFRAME_DIM_ZSTD_LEVEL())

  cfg <- PlatformConfig$new()
  cfg$set('tiledb', 'create', 'dataframe_dim_zstd_level', 999)
  tdco <- TileDBCreateOptions$new(cfg)
  expect_equal(tdco$dataframe_dim_zstd_level(), 999)
})

# Now verify each config value in TileDBCreateOptions.

test_that("TileDBCreateOptions dataframe_dim_zstd_level", {
  cfg <- PlatformConfig$new()
  tdco <- TileDBCreateOptions$new(cfg)
  expect_equal(tdco$dataframe_dim_zstd_level(), DEFAULT_DATAFRAME_DIM_ZSTD_LEVEL())

  cfg <- PlatformConfig$new()
  cfg$set('tiledb', 'create', 'dataframe_dim_zstd_level', 999)
  tdco <- TileDBCreateOptions$new(cfg)
  expect_equal(tdco$dataframe_dim_zstd_level(), 999)
})

test_that("TileDBCreateOptions sparse_nd_array_dim_zstd_level", {
  cfg <- PlatformConfig$new()
  tdco <- TileDBCreateOptions$new(cfg)
  expect_equal(tdco$sparse_nd_array_dim_zstd_level(), DEFAULT_SPARSE_ND_ARRAY_DIM_ZSTD_LEVEL())

  cfg <- PlatformConfig$new()
  cfg$set('tiledb', 'create', 'sparse_nd_array_dim_zstd_level', 999)
  tdco <- TileDBCreateOptions$new(cfg)
  expect_equal(tdco$sparse_nd_array_dim_zstd_level(), 999)
})

test_that("TileDBCreateOptions write_X_chunked", {
  cfg <- PlatformConfig$new()
  tdco <- TileDBCreateOptions$new(cfg)
  expect_equal(tdco$write_X_chunked(), DEFAULT_WRITE_X_CHUNKED())

  cfg <- PlatformConfig$new()
  cfg$set('tiledb', 'create', 'write_X_chunked', FALSE)
  tdco <- TileDBCreateOptions$new(cfg)
  expect_equal(tdco$write_X_chunked(), FALSE)

  cfg <- PlatformConfig$new()
  cfg$set('tiledb', 'create', 'write_X_chunked', TRUE)
  tdco <- TileDBCreateOptions$new(cfg)
  expect_equal(tdco$write_X_chunked(), TRUE)
})

test_that("TileDBCreateOptions goal_chunk_nnz", {
  cfg <- PlatformConfig$new()
  tdco <- TileDBCreateOptions$new(cfg)
  expect_equal(tdco$goal_chunk_nnz(), DEFAULT_GOAL_CHUNK_NNZ())

  cfg <- PlatformConfig$new()
  cfg$set('tiledb', 'create', 'goal_chunk_nnz', 999)
  tdco <- TileDBCreateOptions$new(cfg)
  expect_equal(tdco$goal_chunk_nnz(), 999)
})

test_that("TileDBCreateOptions cell_tile_orders", {
  cfg <- PlatformConfig$new()
  tdco <- TileDBCreateOptions$new(cfg)
  expect_equal(tdco$cell_tile_orders(), c(cell_order = DEFAULT_CELL_ORDER(), tile_order = DEFAULT_TILE_ORDER()))

  cfg <- PlatformConfig$new()
  cfg$set('tiledb', 'create', 'cell_order', 'foo')
  tdco <- TileDBCreateOptions$new(cfg)
  expect_equal(tdco$cell_tile_orders(), c(cell_order = 'foo', tile_order = NULL))

  cfg <- PlatformConfig$new()
  cfg$set('tiledb', 'create', 'tile_order', 'bar')
  tdco <- TileDBCreateOptions$new(cfg)
  expect_equal(tdco$cell_tile_orders(), c(cell_order = NULL, tile_order = 'bar'))

  cfg <- PlatformConfig$new()
  cfg$set('tiledb', 'create', 'cell_order', 'foo')
  cfg$set('tiledb', 'create', 'tile_order', 'bar')
  tdco <- TileDBCreateOptions$new(cfg)
  expect_equal(tdco$cell_tile_orders(), c(cell_order = 'foo', tile_order = 'bar'))
})

test_that("TileDBCreateOptions dim_tile", {
  cfg <- PlatformConfig$new()
  tdco <- TileDBCreateOptions$new(cfg)
  expect_equal(tdco$dim_tile("soma_dim_0"), DEFAULT_TILE_EXTENT())

  cfg <- PlatformConfig$new()
  cfg$set('tiledb', 'create', 'dims', list(soma_dim_0 = list(tile = 999)))
  tdco <- TileDBCreateOptions$new(cfg)
  expect_equal(tdco$dim_tile("soma_dim_0"), 999)
})

# TODO: our ConfigList class currently only accepts scalar values.
# But to port this logic from Python we'll need it to accept list values
# as well.
#
# test_that("TileDBCreateOptions offsets_filters", {
#   cfg <- PlatformConfig$new()
#   tdco <- TileDBCreateOptions$new(cfg)
#   expect_equal(tdco$offsets_filters(), DEFAULT_OFFSETS_FILTERS())
# 
#   cfg <- PlatformConfig$new()
#   cfg$set('tiledb', 'create', 'offsets_filters',
#     list(
#       "RLE",
#       list(name = "ZSTD", COMPRESSION_LEVEL = 9)
#     )
#   )
#   tdco <- TileDBCreateOptions$new(cfg)
#   expect_equal(tdco$offsets_filters(), TBD)
# })

# TODO:
# validity_filters
# dim_filters
# attr_filters

# ================================================================
# TODO: test soma dataframe et al. here, or in their own class-test files
