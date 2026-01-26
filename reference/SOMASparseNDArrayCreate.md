# Create a SOMA Sparse ND Array

Factory function to create a [SOMA sparse ND
array](SOMASparseNDArray.md) for writing (lifecycle: maturing).

## Usage

``` r
SOMASparseNDArrayCreate(
  uri,
  type,
  shape,
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

- type:

  An [Arrow
  type](https://arrow.apache.org/docs/r/reference/data-type.html)
  defining the type of each element in the array.

- shape:

  A vector of integers defining the shape of the array.

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

A new [SOMA sparse ND array](SOMASparseNDArray.md) stored at `uri`
opened for writing.

## Examples

``` r
uri <- withr::local_tempfile(pattern = "soma-sparse-array")
mat <- Matrix::rsparsematrix(100L, 100L, 0.7, repr = "T")
mat[1:3, 1:5]
#> 3 x 5 sparse Matrix of class "dgTMatrix"
#>                                   
#> [1,] -0.073 -2.80 -1.20  .    0.26
#> [2,] -0.560 -0.22 -0.74  0.33 0.28
#> [3,] -0.110 -0.11  .    -0.92 1.50

(arr <- SOMASparseNDArrayCreate(uri, arrow::float64(), shape = dim(mat)))
#> <SOMASparseNDArray>
#>   uri: /tmp/RtmpbAgXbM/soma-sparse-array284656532d03
#>   dimensions: soma_dim_0, soma_dim_1 
#>   attributes: soma_data 
arr$write(mat)
arr$close()

(arr <- SOMASparseNDArrayOpen(uri))
#> <SOMASparseNDArray>
#>   uri: /tmp/RtmpbAgXbM/soma-sparse-array284656532d03
#>   dimensions: soma_dim_0, soma_dim_1 
#>   attributes: soma_data 
m2 <- arr$read()$sparse_matrix()$concat()
m2[1:3, 1:5]
#> 3 x 5 sparse Matrix of class "dgTMatrix"
#>                                   
#> [1,] -0.073 -2.80 -1.20  .    0.26
#> [2,] -0.560 -0.22 -0.74  0.33 0.28
#> [3,] -0.110 -0.11  .    -0.92 1.50
```
