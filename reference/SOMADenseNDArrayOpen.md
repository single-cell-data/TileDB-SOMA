# Open a SOMA Dense Nd Array

Factory function to open a [SOMA dense ND array](SOMADenseNDArray.md)
for reading (lifecycle: maturing).

## Usage

``` r
SOMADenseNDArrayOpen(
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

A [SOMA dense ND array](SOMADenseNDArray.md) stored at `uri` opened in
mode `mode`.

## Examples

``` r
uri <- withr::local_tempfile(pattern = "soma-dense-array")
mat <- matrix(stats::rnorm(100L ^ 2L), nrow = 100L, ncol = 100L)
mat[1:3, 1:5]
#>            [,1]        [,2]       [,3]      [,4]       [,5]
#> [1,]  1.6560563 -0.32979219 -0.7966739 1.4496873  0.1260614
#> [2,] -0.1465172  0.46206288 -1.6059548 0.2504118 -0.3673699
#> [3,]  0.1286809  0.07297064  0.5875243 0.7592450 -0.5737183

(arr <- SOMADenseNDArrayCreate(uri, arrow::float64(), shape = dim(mat)))
#> <SOMADenseNDArray>
#>   uri: /tmp/RtmpbAgXbM/soma-dense-array28466bde969a
#>   dimensions: soma_dim_0, soma_dim_1 
#>   attributes: soma_data 
arr$write(mat)
arr$close()

(arr <- SOMADenseNDArrayOpen(uri))
#> <SOMADenseNDArray>
#>   uri: /tmp/RtmpbAgXbM/soma-dense-array28466bde969a
#>   dimensions: soma_dim_0, soma_dim_1 
#>   attributes: soma_data 
arr$read_arrow_table()
#> Table
#> 10000 rows x 1 columns
#> $soma_data <double not null>
```
