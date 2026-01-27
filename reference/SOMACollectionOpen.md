# Open a SOMA Collection

Factory function to open a [SOMA collection](SOMACollection.md) for
reading (lifecycle: maturing).

## Usage

``` r
SOMACollectionOpen(
  uri,
  mode = "READ",
  platform_config = NULL,
  tiledbsoma_ctx = NULL,
  tiledb_timestamp = NULL,
  context = NULL
)
```

## Arguments

- uri:

  URI for the TileDB object.

- mode:

  One of “`READ`” or “`WRITE`”.

- platform_config:

  Optional platform configuration.

- tiledbsoma_ctx:

  Optional (DEPRECATED) SOMATileDBContext.

- tiledb_timestamp:

  Optional Datetime (POSIXct) for TileDB timestamp; defaults to the
  current time. If not `NULL`, all members accessed through the
  collection inherit the timestamp.

- context:

  Optional `SOMAContext` object used for TileDB operations. If a context
  is not provided, then the default context will be used. Call
  `set_default_context` once before other SOMA operations to configure
  the default context.

## Value

A [SOMA collection](SOMACollection.md) stored at `uri` opened in mode
`mode`.

## Examples

``` r
uri <- withr::local_tempfile(pattern = "soma-collection")

(col <- SOMACollectionCreate(uri))
#> <SOMACollection>
#>   uri: /tmp/Rtmpfr8mYm/soma-collection2a7e1ef348d5
col$add_new_sparse_ndarray("sparse", arrow::float64(), shape = c(100L, 100L))
#> <SOMASparseNDArray>
#>   uri: /tmp/Rtmpfr8mYm/soma-collection2a7e1ef348d5/sparse
#>   dimensions: soma_dim_0, soma_dim_1 
#>   attributes: soma_data 
col$close()

(col <- SOMACollectionOpen(uri))
#> <SOMACollection>
#>   uri: /tmp/Rtmpfr8mYm/soma-collection2a7e1ef348d5
col$names()
#> [1] "sparse"
```
