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
  expect_equal(length(tdco$offsets_filters()), length(DEFAULT_OFFSETS_FILTERS()))

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
  expect_equal(length(tdco$validity_filters()), length(DEFAULT_VALIDITY_FILTERS()))

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
  expect_equal(tdco$dim_tile('soma_dim_1'), DEFAULT_TILE_EXTENT())

  expect_error(tdco$dim_filters())
  expect_equal(length(tdco$dim_filters("soma_dim_0")), 0)
  expect_equal(length(tdco$dim_filters("soma_dim_1")), 0)
  expect_equal(length(tdco$dim_filters("soma_dim_2")), 2)

  expect_error(tdco$attr_filters())
  expect_equal(length(tdco$attr_filters("soma_data_a")), 2)
  expect_equal(length(tdco$attr_filters("soma_data_b")), 0)
})

test_that("platform_config is respected", {
  uri <- withr::local_tempdir("soma-dataframe")

  # Set Arrow schema
  asch <- arrow::schema(
    arrow::field("soma_joinid", arrow::int64(), nullable = FALSE),
    arrow::field("i32", arrow::int32(), nullable = FALSE),
    arrow::field("f64", arrow::float64(), nullable = FALSE),
    arrow::field("utf8", arrow::large_utf8(), nullable = FALSE)
  )

  # Set tiledb create options
  cfg <- PlatformConfig$new()
  cfg$set('tiledb', 'create', 'dataframe_dim_zstd_level', 8)
  cfg$set('tiledb', 'create', 'sparse_nd_array_dim_zstd_level', 9)
  cfg$set('tiledb', 'create', 'capacity', 8000)
  cfg$set('tiledb', 'create', 'tile_order', 'COL_MAJOR')
  cfg$set('tiledb', 'create', 'cell_order', 'UNORDERED')
  cfg$set('tiledb', 'create', 'offsets_filters', list("RLE"))
  cfg$set('tiledb', 'create', 'validity_filters', list("RLE", "NONE"))
  cfg$set('tiledb', 'create', 'dims', list(
    soma_joinid = list(
      filters = list("RLE", list(name="ZSTD", COMPRESSION_LEVEL=8), "NONE")
      # TODO: test setting/checking tile extent, once shapes/domain-maxes are made programmable.
      # At present we get:
      #
      #   Error: Tile extent check failed; domain max expanded to multiple of tile extent exceeds
      #   max value representable by domain type
      #
      # tile = 999
    )
  ))
  cfg$set('tiledb', 'create', 'attrs', list(
    i32 = list(
      filters = list("RLE", list(name="ZSTD", COMPRESSION_LEVEL=9))
    ),
    f64 = list(
      filters = list()
    )
  ))

  # Create the SOMADataFrame
  sdf <- SOMADataFrameCreate(uri=uri, schema=asch, index_column_names=c("soma_joinid"), platform_config = cfg)

  # Read back and check the array schema against the tiledb create options
  arr <- tiledb::tiledb_array(uri)
  tsch <- tiledb::schema(arr)

  # TODO: zstd levels x 2

  expect_equal(tiledb::capacity(tsch), 8000)
  expect_equal(tiledb::tile_order(tsch), "COL_MAJOR")
  expect_equal(tiledb::cell_order(tsch), "UNORDERED")

  offsets_filters <- tiledb::filter_list(tsch)$offsets
  expect_equal(tiledb::nfilters(offsets_filters), 1)
  o1 <- offsets_filters[0] # C++ indexing here
  expect_equal(tiledb::tiledb_filter_type(o1), "RLE")

  validity_filters <- tiledb::filter_list(tsch)$validity
  expect_equal(tiledb::nfilters(validity_filters), 2)
  v1 <- validity_filters[0] # C++ indexing here
  v2 <- validity_filters[1] # C++ indexing here
  expect_equal(tiledb::tiledb_filter_type(v1), "RLE")
  expect_equal(tiledb::tiledb_filter_type(v2), "NONE")

  dom <- tiledb::domain(tsch)
  expect_equal(tiledb::tiledb_ndim(dom), 1)
  dim <- tiledb::dimensions(dom)[[1]]
  expect_equal(tiledb::name(dim), "soma_joinid")
  # TODO: As noted above, check this when we are able to.
  # expect_equal(tiledb::tile(dim), 999)
  dim_filters <- tiledb::filter_list(dim)
  expect_equal(tiledb::nfilters(dim_filters), 3)
  d1 <- dim_filters[0] # C++ indexing here
  d2 <- dim_filters[1] # C++ indexing here
  d3 <- dim_filters[2] # C++ indexing here
  expect_equal(tiledb::tiledb_filter_type(d1), "RLE")
  expect_equal(tiledb::tiledb_filter_type(d2), "ZSTD")
  expect_equal(tiledb::tiledb_filter_type(d3), "NONE")
  expect_equal(tiledb::tiledb_filter_get_option(d2, "COMPRESSION_LEVEL"), 8)

  expect_equal(length(tiledb::attrs(tsch)), 3)
  i32_filters <- tiledb::filter_list(tiledb::attrs(tsch)$i32)
  f64_filters <- tiledb::filter_list(tiledb::attrs(tsch)$f64)
  expect_equal(tiledb::nfilters(i32_filters), 2)
  expect_equal(tiledb::nfilters(f64_filters), 0)

  i1 <- i32_filters[0] # C++ indexing here
  i2 <- i32_filters[1] # C++ indexing here
  expect_equal(tiledb::tiledb_filter_type(i1), "RLE")
  expect_equal(tiledb::tiledb_filter_type(i2), "ZSTD")
  expect_equal(tiledb::tiledb_filter_get_option(i2, "COMPRESSION_LEVEL"), 9)
})
