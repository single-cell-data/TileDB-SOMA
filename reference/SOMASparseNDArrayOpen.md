# Open a SOMA Sparse ND Array

Factory function to open a [SOMA sparse ND array](SOMASparseNDArray.md)
for reading (lifecycle: maturing).

## Usage

``` r
SOMASparseNDArrayOpen(
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
  current time.

- context:

  Optional `SOMAContext` object used for TileDB operations. If a context
  is not provided, then the default context will be used. Call
  `set_default_context` once before other SOMA operations to configure
  the default context.

## Value

A [SOMA sparse ND array](SOMASparseNDArray.md) stored at `uri` opened in
mode `mode`.

## Examples

``` r
uri <- withr::local_tempfile(pattern = "soma-sparse-array")
mat <- Matrix::rsparsematrix(100L, 100L, 0.7, repr = "T")
mat[1:3, 1:5]
#> 3 x 5 sparse Matrix of class "dgTMatrix"
#>                                
#> [1,]  .    -0.42 -0.28 .   0.68
#> [2,]  .    -0.64  1.10 .   .   
#> [3,] -0.47 -0.35 -0.91 1.2 .   

(arr <- SOMASparseNDArrayCreate(uri, arrow::float64(), shape = dim(mat)))
#> <SOMASparseNDArray>
#>   uri: /tmp/RtmpbAgXbM/soma-sparse-array284633198c01
#>   dimensions: soma_dim_0, soma_dim_1 
#>   attributes: soma_data 
arr$write(mat)
arr$close()

(arr <- SOMASparseNDArrayOpen(uri))
#> <SOMASparseNDArray>
#>   uri: /tmp/RtmpbAgXbM/soma-sparse-array284633198c01
#>   dimensions: soma_dim_0, soma_dim_1 
#>   attributes: soma_data 
m2 <- arr$read()$sparse_matrix()$concat()
m2[1:3, 1:5]
#> 3 x 5 sparse Matrix of class "dgTMatrix"
#>                                
#> [1,]  .    -0.42 -0.28 .   0.68
#> [2,]  .    -0.64  1.10 .   .   
#> [3,] -0.47 -0.35 -0.91 1.2 .   
```
