# TileDBCreateOptions

Provides strongly-typed access and default values for `platform_config`
options stored under the “tiledb” \\\rightarrow\\ “create” mapping keys.

Intended for internal use only.

## Super class

[`tiledbsoma::MappingBase`](MappingBase.md) -\> `TileDBCreateOptions`

## Methods

### Public methods

- [`TileDBCreateOptions$new()`](#method-TileDBCreateOptions-new)

- [`TileDBCreateOptions$cell_tile_orders()`](#method-TileDBCreateOptions-cell_tile_orders)

- [`TileDBCreateOptions$dim_tile()`](#method-TileDBCreateOptions-dim_tile)

- [`TileDBCreateOptions$capacity()`](#method-TileDBCreateOptions-capacity)

- [`TileDBCreateOptions$allows_duplicates()`](#method-TileDBCreateOptions-allows_duplicates)

- [`TileDBCreateOptions$dataframe_dim_zstd_level()`](#method-TileDBCreateOptions-dataframe_dim_zstd_level)

- [`TileDBCreateOptions$sparse_nd_array_dim_zstd_level()`](#method-TileDBCreateOptions-sparse_nd_array_dim_zstd_level)

- [`TileDBCreateOptions$dense_nd_array_dim_zstd_level()`](#method-TileDBCreateOptions-dense_nd_array_dim_zstd_level)

- [`TileDBCreateOptions$offsets_filters()`](#method-TileDBCreateOptions-offsets_filters)

- [`TileDBCreateOptions$validity_filters()`](#method-TileDBCreateOptions-validity_filters)

- [`TileDBCreateOptions$dim_filters()`](#method-TileDBCreateOptions-dim_filters)

- [`TileDBCreateOptions$attr_filters()`](#method-TileDBCreateOptions-attr_filters)

- [`TileDBCreateOptions$write_X_chunked()`](#method-TileDBCreateOptions-write_X_chunked)

- [`TileDBCreateOptions$goal_chunk_nnz()`](#method-TileDBCreateOptions-goal_chunk_nnz)

- [`TileDBCreateOptions$to_list()`](#method-TileDBCreateOptions-to_list)

- [`TileDBCreateOptions$clone()`](#method-TileDBCreateOptions-clone)

Inherited methods

- [`tiledbsoma::MappingBase$get()`](MappingBase.html#method-get)
- [`tiledbsoma::MappingBase$items()`](MappingBase.html#method-items)
- [`tiledbsoma::MappingBase$keys()`](MappingBase.html#method-keys)
- [`tiledbsoma::MappingBase$length()`](MappingBase.html#method-length)
- [`tiledbsoma::MappingBase$print()`](MappingBase.html#method-print)
- [`tiledbsoma::MappingBase$remove()`](MappingBase.html#method-remove)
- [`tiledbsoma::MappingBase$set()`](MappingBase.html#method-set)
- [`tiledbsoma::MappingBase$setv()`](MappingBase.html#method-setv)
- [`tiledbsoma::MappingBase$update()`](MappingBase.html#method-update)
- [`tiledbsoma::MappingBase$values()`](MappingBase.html#method-values)

------------------------------------------------------------------------

### Method `new()`

Create a `TileDBCreateOptions` object

#### Usage

    TileDBCreateOptions$new(platform_config = NULL)

#### Arguments

- `platform_config`:

  A [platform configuration](PlatformConfig.md) object

------------------------------------------------------------------------

### Method `cell_tile_orders()`

Returns the cell and tile orders that should be used. If neither
`cell_order` nor `tile_order` is present, only in this case will we use
the default values provided.

#### Usage

    TileDBCreateOptions$cell_tile_orders()

#### Returns

A two-length character vector with names of “`cell_order`” and
“`tile_order`”

------------------------------------------------------------------------

### Method `dim_tile()`

#### Usage

    TileDBCreateOptions$dim_tile(dim_name, default = 2048)

#### Arguments

- `dim_name`:

  Name of dimension to get tiling for

- `default`:

  Default tiling if `dim_name` is not set

#### Returns

int

#### Examples

    cfg <- PlatformConfig$new()
    cfg$set(
      platform = 'tiledb',
      param = 'create',
      key = 'dims',
      value = list(soma_dim_0 = list(tile = 999))
    )
    (tdco <- TileDBCreateOptions$new(cfg))
    tdco$dim_tile("soma_dim_0")
    tdco$dim_tile("soma_dim_1")

------------------------------------------------------------------------

### Method `capacity()`

#### Usage

    TileDBCreateOptions$capacity()

#### Returns

int

------------------------------------------------------------------------

### Method `allows_duplicates()`

#### Usage

    TileDBCreateOptions$allows_duplicates()

#### Returns

bool

------------------------------------------------------------------------

### Method `dataframe_dim_zstd_level()`

#### Usage

    TileDBCreateOptions$dataframe_dim_zstd_level()

#### Returns

int

------------------------------------------------------------------------

### Method `sparse_nd_array_dim_zstd_level()`

#### Usage

    TileDBCreateOptions$sparse_nd_array_dim_zstd_level()

#### Returns

int

------------------------------------------------------------------------

### Method `dense_nd_array_dim_zstd_level()`

#### Usage

    TileDBCreateOptions$dense_nd_array_dim_zstd_level()

#### Returns

int

------------------------------------------------------------------------

### Method `offsets_filters()`

#### Usage

    TileDBCreateOptions$offsets_filters(default = list())

#### Arguments

- `default`:

  Default offset filters to use if not currently set

#### Returns

A list of
[`tiledb_filter`](https://tiledb-inc.github.io/TileDB-R/reference/tiledb_filter-class.html)
objects

------------------------------------------------------------------------

### Method `validity_filters()`

#### Usage

    TileDBCreateOptions$validity_filters(default = list())

#### Arguments

- `default`:

  Default validity filters to use if not currently set

#### Returns

A list of
[`tiledb_filter`](https://tiledb-inc.github.io/TileDB-R/reference/tiledb_filter-class.html)
objects

------------------------------------------------------------------------

### Method `dim_filters()`

#### Usage

    TileDBCreateOptions$dim_filters(dim_name, default = list())

#### Arguments

- `dim_name`:

  Name of dimension to get filters for

- `default`:

  Default filters to use for if not currently set

#### Returns

A list of
[`tiledb_filter`](https://tiledb-inc.github.io/TileDB-R/reference/tiledb_filter-class.html)
objects

#### Examples

    filters <- list(
      soma_dim_0 = list(tile = 100, filters = list("RLE")),
      soma_dim_1 = list(tile = 200, filters = list("RLE", list(name = "ZSTD", COMPRESSION_LEVEL = 9)))
    )
    cfg <- PlatformConfig$new()
    cfg$set(platform = 'tiledb', param = 'create', key = 'dims', value = filters)
    (tdco <- TileDBCreateOptions$new(cfg))
    tdco$dim_filters("soma_dim_0")
    tdco$dim_filters("non-existant")

------------------------------------------------------------------------

### Method `attr_filters()`

#### Usage

    TileDBCreateOptions$attr_filters(attr_name, default = list())

#### Arguments

- `attr_name`:

  Name of attribute

- `default`:

  Default filters to use if not currently set

#### Returns

A list of
[`tiledb_filter`](https://tiledb-inc.github.io/TileDB-R/reference/tiledb_filter-class.html)
objects

#### Examples

    filters <- list(
      soma_data_a = list(filters = list("RLE")),
      soma_data_b = list(filters = list("RLE", list(name = "ZSTD", COMPRESSION_LEVEL = 9)))
    )
    cfg <- PlatformConfig$new()
    cfg$set(platform = 'tiledb', param = 'create', key = 'attrs', value = filters)
    (tdco <- TileDBCreateOptions$new(cfg))
    tdco$attr_filters("soma_data_b")
    tdco$attr_filters("non-existant")

------------------------------------------------------------------------

### Method `write_X_chunked()`

#### Usage

    TileDBCreateOptions$write_X_chunked()

#### Returns

bool

------------------------------------------------------------------------

### Method `goal_chunk_nnz()`

#### Usage

    TileDBCreateOptions$goal_chunk_nnz()

#### Returns

int

------------------------------------------------------------------------

### Method `to_list()`

...

#### Usage

    TileDBCreateOptions$to_list(build_filters = TRUE)

#### Arguments

- `build_filters`:

  Build filters into
  [`tiledb_filter`](https://tiledb-inc.github.io/TileDB-R/reference/tiledb_filter-class.html)
  objects. If set to `FALSE`, JSON strings are created instead of filter
  objects.

#### Returns

The 'create options' as a list

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    TileDBCreateOptions$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
(cfg <- PlatformConfig$new())
#> <PlatformConfig>
(tdco <- TileDBCreateOptions$new(cfg))
#> <TileDBCreateOptions>
#>   tile_order: ROW_MAJOR
#>   cell_order: ROW_MAJOR
#>   capacity: 1e+05
#>   allows_duplicates: FALSE
#>   dataframe_dim_zstd_level: 3
#>   sparse_nd_array_dim_zstd_level: 3
#>   dense_nd_array_dim_zstd_level: 3
#>   offsets_filters: list("DOUBLE_DELTA", "BIT_WIDTH_REDUCTION", "ZSTD")
#>   write_X_chunked: TRUE
#>   goal_chunk_nnz: 2e+08
tdco$cell_tile_orders()
#>  cell_order  tile_order 
#> "ROW_MAJOR" "ROW_MAJOR" 
tdco$to_list()
#> $tile_order
#> [1] "ROW_MAJOR"
#> 
#> $cell_order
#> [1] "ROW_MAJOR"
#> 
#> $capacity
#> [1] 1e+05
#> 
#> $allows_duplicates
#> [1] FALSE
#> 
#> $dataframe_dim_zstd_level
#> [1] 3
#> 
#> $sparse_nd_array_dim_zstd_level
#> [1] 3
#> 
#> $dense_nd_array_dim_zstd_level
#> [1] 3
#> 
#> $offsets_filters
#> $offsets_filters[[1]]
#> tiledb_filter("DOUBLE_DELTA") 
#> 
#> $offsets_filters[[2]]
#> tiledb_filter_set_option(tiledb_filter("BIT_WIDTH_REDUCTION"),"BIT_WIDTH_MAX_WINDOW",256) 
#> 
#> $offsets_filters[[3]]
#> tiledb_filter_set_option(tiledb_filter("ZSTD"),"COMPRESSION_LEVEL",-1) 
#> 
#> 
#> $write_X_chunked
#> [1] TRUE
#> 
#> $goal_chunk_nnz
#> [1] 2e+08
#> 
#> $validity_filters
#> list()
#> 
tdco$to_list(build_filters = FALSE)
#> $tile_order
#> [1] "ROW_MAJOR"
#> 
#> $cell_order
#> [1] "ROW_MAJOR"
#> 
#> $capacity
#> [1] 1e+05
#> 
#> $allows_duplicates
#> [1] FALSE
#> 
#> $dataframe_dim_zstd_level
#> [1] 3
#> 
#> $sparse_nd_array_dim_zstd_level
#> [1] 3
#> 
#> $dense_nd_array_dim_zstd_level
#> [1] 3
#> 
#> $offsets_filters
#> [1] "[{  \"name\": \"DOUBLE_DELTA\" }, {  \"name\": \"BIT_WIDTH_REDUCTION\" }, {  \"name\": \"ZSTD\" }]"
#> 
#> $write_X_chunked
#> [1] TRUE
#> 
#> $goal_chunk_nnz
#> [1] 2e+08
#> 
#> $validity_filters
#> [1] ""
#> 
#> $dims
#> [1] "{ }"
#> 
#> $attrs
#> [1] "{ }"
#> 


## ------------------------------------------------
## Method `TileDBCreateOptions$dim_tile`
## ------------------------------------------------

cfg <- PlatformConfig$new()
cfg$set(
  platform = 'tiledb',
  param = 'create',
  key = 'dims',
  value = list(soma_dim_0 = list(tile = 999))
)
(tdco <- TileDBCreateOptions$new(cfg))
#> <TileDBCreateOptions>
#>   dims: list(soma_dim_0 = list(tile = 999))
#>   tile_order: ROW_MAJOR
#>   cell_order: ROW_MAJOR
#>   capacity: 1e+05
#>   allows_duplicates: FALSE
#>   dataframe_dim_zstd_level: 3
#>   sparse_nd_array_dim_zstd_level: 3
#>   dense_nd_array_dim_zstd_level: 3
#>   offsets_filters: list("DOUBLE_DELTA", "BIT_WIDTH_REDUCTION", "ZSTD")
#>   write_X_chunked: TRUE
#>   goal_chunk_nnz: 2e+08
tdco$dim_tile("soma_dim_0")
#> [1] 999
tdco$dim_tile("soma_dim_1")
#> [1] 2048


## ------------------------------------------------
## Method `TileDBCreateOptions$dim_filters`
## ------------------------------------------------

filters <- list(
  soma_dim_0 = list(tile = 100, filters = list("RLE")),
  soma_dim_1 = list(tile = 200, filters = list("RLE", list(name = "ZSTD", COMPRESSION_LEVEL = 9)))
)
cfg <- PlatformConfig$new()
cfg$set(platform = 'tiledb', param = 'create', key = 'dims', value = filters)
(tdco <- TileDBCreateOptions$new(cfg))
#> <TileDBCreateOptions>
#>   dims: list(soma_dim_0 = list(tile = 100, filters = list("RLE")), soma_dim_1 = list(tile = 200, filters = list("RLE", list(name = "ZSTD", COMPRESSION_LEVEL = 9))))
#>   tile_order: ROW_MAJOR
#>   cell_order: ROW_MAJOR
#>   capacity: 1e+05
#>   allows_duplicates: FALSE
#>   dataframe_dim_zstd_level: 3
#>   sparse_nd_array_dim_zstd_level: 3
#>   dense_nd_array_dim_zstd_level: 3
#>   offsets_filters: list("DOUBLE_DELTA", "BIT_WIDTH_REDUCTION", "ZSTD")
#>   write_X_chunked: TRUE
#>   goal_chunk_nnz: 2e+08
tdco$dim_filters("soma_dim_0")
#> [[1]]
#> tiledb_filter_set_option(tiledb_filter("RLE"),"COMPRESSION_LEVEL",-1) 
#> 
tdco$dim_filters("non-existant")
#> list()


## ------------------------------------------------
## Method `TileDBCreateOptions$attr_filters`
## ------------------------------------------------

filters <- list(
  soma_data_a = list(filters = list("RLE")),
  soma_data_b = list(filters = list("RLE", list(name = "ZSTD", COMPRESSION_LEVEL = 9)))
)
cfg <- PlatformConfig$new()
cfg$set(platform = 'tiledb', param = 'create', key = 'attrs', value = filters)
(tdco <- TileDBCreateOptions$new(cfg))
#> <TileDBCreateOptions>
#>   attrs: list(soma_data_a = list(filters = list("RLE")), soma_data_b = list(filters = list("RLE", list(name = "ZSTD", COMPRESSION_LEVEL = 9))))
#>   tile_order: ROW_MAJOR
#>   cell_order: ROW_MAJOR
#>   capacity: 1e+05
#>   allows_duplicates: FALSE
#>   dataframe_dim_zstd_level: 3
#>   sparse_nd_array_dim_zstd_level: 3
#>   dense_nd_array_dim_zstd_level: 3
#>   offsets_filters: list("DOUBLE_DELTA", "BIT_WIDTH_REDUCTION", "ZSTD")
#>   write_X_chunked: TRUE
#>   goal_chunk_nnz: 2e+08
tdco$attr_filters("soma_data_b")
#> [[1]]
#> tiledb_filter_set_option(tiledb_filter("RLE"),"COMPRESSION_LEVEL",-1) 
#> 
#> [[2]]
#> tiledb_filter_set_option(tiledb_filter("ZSTD"),"COMPRESSION_LEVEL",9) 
#> 
tdco$attr_filters("non-existant")
#> list()
```
