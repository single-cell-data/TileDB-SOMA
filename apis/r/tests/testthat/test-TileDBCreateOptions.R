# Verify that the TileDBCreateOptions takes a NULL or a PlatformConfig
test_that("TileDBCreateOptions construction", {
  expect_error(TileDBCreateOptions$new(c(foo = "bar")))

  tdco <- TileDBCreateOptions$new(NULL)

  cfg <- PlatformConfig$new()
  tdco <- TileDBCreateOptions$new(cfg)
})

# Verify access of parameters across a range of nestedness of keys in the PlatformConfig.
test_that("TileDBCreateOptions access from PlatformConfig", {
  cfg <- PlatformConfig$new()
  tdco <- TileDBCreateOptions$new(cfg)
  expect_equal(tdco$dataframe_dim_zstd_level(), .CREATE_DEFAULTS$dataframe_dim_zstd_level)

  cfg <- PlatformConfig$new()
  cfg$set('not_tiledb', 'not_create', 'not_dataframe_dim_zstd_level', 999)
  tdco <- TileDBCreateOptions$new(cfg)
  expect_equal(tdco$dataframe_dim_zstd_level(), .CREATE_DEFAULTS$dataframe_dim_zstd_level)

  cfg <- PlatformConfig$new()
  cfg$set('tiledb', 'not_create', 'not_dataframe_dim_zstd_level', 999)
  tdco <- TileDBCreateOptions$new(cfg)
  expect_equal(tdco$dataframe_dim_zstd_level(), .CREATE_DEFAULTS$dataframe_dim_zstd_level)

  cfg <- PlatformConfig$new()
  cfg$set('not_tiledb', 'create', 'not_dataframe_dim_zstd_level', 999)
  tdco <- TileDBCreateOptions$new(cfg)
  expect_equal(tdco$dataframe_dim_zstd_level(), .CREATE_DEFAULTS$dataframe_dim_zstd_level)

  cfg <- PlatformConfig$new()
  cfg$set('tiledb', 'create', 'not_dataframe_dim_zstd_level', 999)
  tdco <- TileDBCreateOptions$new(cfg)
  expect_equal(tdco$dataframe_dim_zstd_level(), .CREATE_DEFAULTS$dataframe_dim_zstd_level)

  cfg <- PlatformConfig$new()
  cfg$set('tiledb', 'create', 'dataframe_dim_zstd_level', 999)
  tdco <- TileDBCreateOptions$new(cfg)
  expect_equal(tdco$dataframe_dim_zstd_level(), 999)
})

# Now verify each config value in TileDBCreateOptions.

test_that("TileDBCreateOptions dataframe_dim_zstd_level", {
  cfg <- PlatformConfig$new()
  tdco <- TileDBCreateOptions$new(cfg)
  expect_equal(tdco$dataframe_dim_zstd_level(), .CREATE_DEFAULTS$dataframe_dim_zstd_level)

  cfg <- PlatformConfig$new()
  cfg$set('tiledb', 'create', 'dataframe_dim_zstd_level', 999)
  tdco <- TileDBCreateOptions$new(cfg)
  expect_equal(tdco$dataframe_dim_zstd_level(), 999)
})

test_that("TileDBCreateOptions sparse_nd_array_dim_zstd_level", {
  cfg <- PlatformConfig$new()
  tdco <- TileDBCreateOptions$new(cfg)
  expect_equal(
    tdco$sparse_nd_array_dim_zstd_level(),
    .CREATE_DEFAULTS$sparse_nd_array_dim_zstd_level
  )

  cfg <- PlatformConfig$new()
  cfg$set('tiledb', 'create', 'sparse_nd_array_dim_zstd_level', 999)
  tdco <- TileDBCreateOptions$new(cfg)
  expect_equal(tdco$sparse_nd_array_dim_zstd_level(), 999)
})

test_that("TileDBCreateOptions write_X_chunked", {
  cfg <- PlatformConfig$new()
  tdco <- TileDBCreateOptions$new(cfg)
  expect_equal(tdco$write_X_chunked(), .CREATE_DEFAULTS$write_X_chunked)

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
  expect_equal(tdco$goal_chunk_nnz(), .CREATE_DEFAULTS$goal_chunk_nnz)

  cfg <- PlatformConfig$new()
  cfg$set('tiledb', 'create', 'goal_chunk_nnz', 999)
  tdco <- TileDBCreateOptions$new(cfg)
  expect_equal(tdco$goal_chunk_nnz(), 999)
})

test_that("TileDBCreateOptions cell_tile_orders", {
  cfg <- PlatformConfig$new()
  tdco <- TileDBCreateOptions$new(cfg)
  expect_equal(
    tdco$cell_tile_orders(),
    c(cell_order = .CREATE_DEFAULTS$cell_order, tile_order = .CREATE_DEFAULTS$tile_order)
  )

  cfg <- PlatformConfig$new()
  cfg$set('tiledb', 'create', 'cell_order', 'foo')
  tdco <- TileDBCreateOptions$new(cfg)
  expect_equal(
    tdco$cell_tile_orders(),
    c(cell_order = 'foo', tile_order = .CREATE_DEFAULTS$tile_order)
  )
  # expect_equal(tdco$cell_tile_orders(), c(cell_order = 'foo', tile_order = NULL))

  cfg <- PlatformConfig$new()
  cfg$set('tiledb', 'create', 'tile_order', 'bar')
  tdco <- TileDBCreateOptions$new(cfg)
  expect_equal(
    tdco$cell_tile_orders(),
    c(cell_order = .CREATE_DEFAULTS$cell_order, tile_order = 'bar')
  )
  # expect_equal(tdco$cell_tile_orders(), c(cell_order = NULL, tile_order = 'bar'))

  cfg <- PlatformConfig$new()
  cfg$set('tiledb', 'create', 'cell_order', 'foo')
  cfg$set('tiledb', 'create', 'tile_order', 'bar')
  tdco <- TileDBCreateOptions$new(cfg)
  expect_equal(tdco$cell_tile_orders(), c(cell_order = 'foo', tile_order = 'bar'))
})

test_that("TileDBCreateOptions dim_tile", {
  cfg <- PlatformConfig$new()
  tdco <- TileDBCreateOptions$new(cfg)
  expect_equal(tdco$dim_tile("soma_dim_0"), 2048)

  cfg <- PlatformConfig$new()
  cfg$set('tiledb', 'create', 'dims', list(soma_dim_0 = list(tile = 999)))
  tdco <- TileDBCreateOptions$new(cfg)
  expect_equal(tdco$dim_tile("soma_dim_0"), 999)

  expect_error(tdco$dim_tile())
})

test_that("TileDBCreateOptions dim_filters", {
  cfg <- PlatformConfig$new()
  tdco <- TileDBCreateOptions$new(cfg)
  expect_error(tdco$dim_filters())
  expect_no_condition(length(tdco$dim_filters("soma_dim_0", default=list("ZSTD"))))

  cfg <- PlatformConfig$new()
  cfg$set('tiledb', 'create', 'dims', list(
    soma_dim_0 = list(filters = list("RLE")),
    soma_dim_1 = list(filters = list("RLE", list(name="ZSTD", COMPRESSION_LEVEL=9)))
  ))
  tdco <- TileDBCreateOptions$new(cfg)
  expect_equal(length(tdco$dim_filters("soma_dim_0", default=list("ZSTD"))), 1)
  expect_equal(length(tdco$dim_filters("soma_dim_1", default=list("ZSTD"))), 2)
})

test_that("TileDBCreateOptions attr_filters", {
  cfg <- PlatformConfig$new()
  tdco <- TileDBCreateOptions$new(cfg)
  expect_error(tdco$attr_filters())
  expect_no_condition(length(tdco$attr_filters("soma_data")))

  cfg <- PlatformConfig$new()
  cfg$set('tiledb', 'create', 'attrs', list(
    soma_data_a = list(filters = list("RLE")),
    soma_data_b = list(filters = list("RLE", list(name="ZSTD", COMPRESSION_LEVEL=9)))
  ))
  tdco <- TileDBCreateOptions$new(cfg)
  expect_error(tdco$attr_filters())
  expect_equal(length(tdco$attr_filters("soma_data_a")), 1)
  expect_equal(length(tdco$attr_filters("soma_data_b")), 2)
  expect_equal(length(tdco$attr_filters("soma_data_c")), 0)
  expect_equal(length(tdco$attr_filters("soma_data_c", list())), 0)
})

test_that("TileDBCreateOptions offsets_filters", {
  cfg <- PlatformConfig$new()
  tdco <- TileDBCreateOptions$new(cfg)
  expect_length(tdco$offsets_filters(), length(.CREATE_DEFAULTS$offsets_filters))

  cfg <- PlatformConfig$new()
  cfg$set('tiledb', 'create', 'offsets_filters',
    list("RLE")
  )
  tdco <- TileDBCreateOptions$new(cfg)
  expect_equal(length(tdco$offsets_filters()), 1)

  cfg <- PlatformConfig$new()
  cfg$set('tiledb', 'create', 'offsets_filters',
    list(
      "RLE",
      list(name = "ZSTD", COMPRESSION_LEVEL = 9)
    )
  )
  tdco <- TileDBCreateOptions$new(cfg)
  expect_equal(length(tdco$offsets_filters()), 2)
})

test_that("TileDBCreateOptions validity_filters", {
  cfg <- PlatformConfig$new()
  tdco <- TileDBCreateOptions$new(cfg)
  expect_length(tdco$validity_filters(), length(.CREATE_DEFAULTS$validity_filters))

  cfg <- PlatformConfig$new()
  cfg$set('tiledb', 'create', 'validity_filters',
    list("RLE")
  )
  tdco <- TileDBCreateOptions$new(cfg)
  expect_equal(length(tdco$validity_filters()), 1)

  cfg <- PlatformConfig$new()
  cfg$set('tiledb', 'create', 'validity_filters',
    list(
      "RLE",
      list(name = "ZSTD", COMPRESSION_LEVEL = 9)
    )
  )
  tdco <- TileDBCreateOptions$new(cfg)
  expect_equal(length(tdco$validity_filters()), 2)
})

test_that("TileDBCreateOptions overrides", {
  cfg <- PlatformConfig$new()
  cfg$set('tiledb', 'create', 'dataframe_dim_zstd_level', 8)
  cfg$set('tiledb', 'create', 'sparse_nd_array_dim_zstd_level', 9)
  cfg$set('tiledb', 'create', 'dims', list(
    soma_dim_0 = list(tile = 6),
    soma_dim_1 = list(filters = list()),
    soma_dim_2 = list(filters = list("RLE", "ZSTD"))
  ))
  cfg$set('tiledb', 'create', 'attrs', list(
    soma_data_a = list(filters = list("RLE", list(name="ZSTD", COMPRESSION_LEVEL=9)))
  ))
  tdco <- TileDBCreateOptions$new(cfg)

  expect_equal(tdco$dataframe_dim_zstd_level(), 8)
  expect_equal(tdco$sparse_nd_array_dim_zstd_level(), 9)

  expect_error(tdco$dim_tile())
  expect_equal(tdco$dim_tile('soma_dim_0'), 6)
  expect_equal(tdco$dim_tile('soma_dim_1'), 2048)

  expect_error(tdco$dim_filters())
  expect_equal(length(tdco$dim_filters("soma_dim_0")), 0)
  expect_equal(length(tdco$dim_filters("soma_dim_1")), 0)
  expect_equal(length(tdco$dim_filters("soma_dim_2")), 2)

  expect_error(tdco$attr_filters())
  expect_equal(length(tdco$attr_filters("soma_data_a")), 2)
  expect_equal(length(tdco$attr_filters("soma_data_b")), 0)
})
