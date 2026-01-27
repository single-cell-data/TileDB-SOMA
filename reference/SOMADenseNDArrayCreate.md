# Create a SOMA Dense ND Array

Factory function to create a [SOMA dense ND array](SOMADenseNDArray.md)
for writing (lifecycle: maturing).

## Usage

``` r
SOMADenseNDArrayCreate(
  uri,
  type,
  shape,
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

A new [SOMA dense ND array](SOMADenseNDArray.md) stored at `uri` opened
for writing.

## Examples

``` r
uri <- withr::local_tempfile(pattern = "soma-dense-array")
mat <- matrix(stats::rnorm(100L ^ 2L), nrow = 100L, ncol = 100L)
mat[1:3, 1:5]
#>             [,1]       [,2]       [,3]      [,4]      [,5]
#> [1,]  0.03052902 -1.2368860 -0.6614749 0.5568221 0.1034234
#> [2,] -0.04787012 -0.5746649 -1.2197908 0.5519798 0.1276377
#> [3,]  1.59039312  0.1616978 -0.4527197 0.6139520 2.0793888

(arr <- SOMADenseNDArrayCreate(uri, arrow::float64(), shape = dim(mat)))
#> <SOMADenseNDArray>
#>   uri: /tmp/Rtmpfr8mYm/soma-dense-array2a7e791df1cc
#>   dimensions: soma_dim_0, soma_dim_1 
#>   attributes: soma_data 
arr$write(mat)
arr$close()

(arr <- SOMADenseNDArrayOpen(uri))
#> <SOMADenseNDArray>
#>   uri: /tmp/Rtmpfr8mYm/soma-dense-array2a7e791df1cc
#>   dimensions: soma_dim_0, soma_dim_1 
#>   attributes: soma_data 
arr$read_arrow_table()
#> Table
#> 10000 rows x 1 columns
#> $soma_data <double not null>
```
