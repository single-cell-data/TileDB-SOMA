# Create a SOMA Collection

Factory function to create a [SOMA collection](SOMACollection.md) for
writing (lifecycle: maturing).

## Usage

``` r
SOMACollectionCreate(
  uri,
  ingest_mode = c("write", "resume"),
  platform_config = NULL,
  tiledbsoma_ctx = NULL,
  tiledb_timestamp = NULL,
  context = NULL
)
```

## Arguments

- uri:

  URI for the TileDB object.

- ingest_mode:

  Ingestion mode when creating the TileDB object; choose from:

  - “`write`”: create a new TileDB object and error if it already
    exists.

  - “`resume`”: attempt to create a new TileDB object; if it already
    exists, simply open it for writing.

- platform_config:

  Optional platform configuration.

- tiledbsoma_ctx:

  Optional (DEPRECATED) SOMATileDBContext.

- tiledb_timestamp:

  Optional Datetime (POSIXct) for TileDB timestamp.

- context:

  Optional `SOMAContext` object used for TileDB operations. If a context
  is not provided, then the default context will be used. Call
  `set_default_context` once before other SOMA operations to configure
  the default context.

## Value

A new [SOMA collection](SOMACollection.md) stored at `uri` opened for
writing.

## Examples

``` r
uri <- withr::local_tempfile(pattern = "soma-collection")

(col <- SOMACollectionCreate(uri))
#> <SOMACollection>
#>   uri: /tmp/Rtmpfr8mYm/soma-collection2a7e397b1d88
col$add_new_sparse_ndarray("sparse", arrow::float64(), shape = c(100L, 100L))
#> <SOMASparseNDArray>
#>   uri: /tmp/Rtmpfr8mYm/soma-collection2a7e397b1d88/sparse
#>   dimensions: soma_dim_0, soma_dim_1 
#>   attributes: soma_data 
col$close()

(col <- SOMACollectionOpen(uri))
#> <SOMACollection>
#>   uri: /tmp/Rtmpfr8mYm/soma-collection2a7e397b1d88
col$names()
#> [1] "sparse"
```
