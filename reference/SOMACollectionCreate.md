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
  tiledb_timestamp = NULL
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

  Optional SOMATileDBContext.

- tiledb_timestamp:

  Optional Datetime (POSIXct) for TileDB timestamp.

## Value

A new [SOMA collection](SOMACollection.md) stored at `uri` opened for
writing.

## Examples

``` r
uri <- withr::local_tempfile(pattern = "soma-collection")

(col <- SOMACollectionCreate(uri))
#> <SOMACollection>
#>   uri: /tmp/Rtmpqg2whH/soma-collection2b78fc1eed2
col$add_new_sparse_ndarray("sparse", arrow::float64(), shape = c(100L, 100L))
#> <SOMASparseNDArray>
#>   uri: /tmp/Rtmpqg2whH/soma-collection2b78fc1eed2/sparse
#>   dimensions: soma_dim_0, soma_dim_1 
#>   attributes: soma_data 
col$close()

(col <- SOMACollectionOpen(uri))
#> <SOMACollection>
#>   uri: /tmp/Rtmpqg2whH/soma-collection2b78fc1eed2
col$names()
#> [1] "sparse"
```
