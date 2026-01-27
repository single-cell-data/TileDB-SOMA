# SOMA Dense Nd-Array

`SOMADenseNDArray` is a dense, N-dimensional array of a `primitive`
type, with offset (zero-based) `int64` integer indexing on each
dimension with domain `[0, maxInt64)`. The `SOMADenseNDArray` has a
user-defined schema, which includes:

- `type`: a `primitive` type, expressed as an Arrow type (e.g.,
  [`int64`](https://arrow.apache.org/docs/r/reference/data-type.html),
  [`float32`](https://arrow.apache.org/docs/r/reference/data-type.html),
  etc), indicating the type of data contained within the array.

- `shape`: the shape of the array, i.e., number and length of each
  dimension. This is a soft limit which can be increased using
  `$resize()` up to the `maxshape`.

- `maxshape`: the hard limit up to which `shape` may be increased using
  `$resize()`.

All dimensions must have a positive, non-zero length, and there must be
1 or more dimensions.

The default “fill” value for `SOMADenseNDArray` is the zero or null
value of the array type (e.g.,
[`arrow::float32()`](https://arrow.apache.org/docs/r/reference/data-type.html)
defaults to 0.0).

## Super classes

[`tiledbsoma::SOMAObject`](SOMAObject.md) -\>
[`tiledbsoma::SOMAArrayBase`](SOMAArrayBase.md) -\>
[`tiledbsoma::SOMANDArrayBase`](SOMANDArrayBase.md) -\>
`SOMADenseNDArray`

## Methods

### Public methods

- [`SOMADenseNDArray$read_arrow_table()`](#method-SOMADenseNDArray-read_arrow_table)

- [`SOMADenseNDArray$read_dense_matrix()`](#method-SOMADenseNDArray-read_dense_matrix)

- [`SOMADenseNDArray$write()`](#method-SOMADenseNDArray-write)

- [`SOMADenseNDArray$clone()`](#method-SOMADenseNDArray-clone)

Inherited methods

- [`tiledbsoma::SOMAObject$class()`](SOMAObject.html#method-class)
- [`tiledbsoma::SOMAObject$exists()`](SOMAObject.html#method-exists)
- [`tiledbsoma::SOMAObject$get_metadata()`](SOMAObject.html#method-get_metadata)
- [`tiledbsoma::SOMAObject$initialize()`](SOMAObject.html#method-initialize)
- [`tiledbsoma::SOMAObject$is_open()`](SOMAObject.html#method-is_open)
- [`tiledbsoma::SOMAObject$mode()`](SOMAObject.html#method-mode)
- [`tiledbsoma::SOMAObject$reopen()`](SOMAObject.html#method-reopen)
- [`tiledbsoma::SOMAObject$set_metadata()`](SOMAObject.html#method-set_metadata)
- [`tiledbsoma::SOMAArrayBase$allows_duplicates()`](SOMAArrayBase.html#method-allows_duplicates)
- [`tiledbsoma::SOMAArrayBase$attributes()`](SOMAArrayBase.html#method-attributes)
- [`tiledbsoma::SOMAArrayBase$attrnames()`](SOMAArrayBase.html#method-attrnames)
- [`tiledbsoma::SOMAArrayBase$close()`](SOMAArrayBase.html#method-close)
- [`tiledbsoma::SOMAArrayBase$colnames()`](SOMAArrayBase.html#method-colnames)
- [`tiledbsoma::SOMAArrayBase$dimensions()`](SOMAArrayBase.html#method-dimensions)
- [`tiledbsoma::SOMAArrayBase$dimnames()`](SOMAArrayBase.html#method-dimnames)
- [`tiledbsoma::SOMAArrayBase$index_column_names()`](SOMAArrayBase.html#method-index_column_names)
- [`tiledbsoma::SOMAArrayBase$is_sparse()`](SOMAArrayBase.html#method-is_sparse)
- [`tiledbsoma::SOMAArrayBase$maxshape()`](SOMAArrayBase.html#method-maxshape)
- [`tiledbsoma::SOMAArrayBase$ndim()`](SOMAArrayBase.html#method-ndim)
- [`tiledbsoma::SOMAArrayBase$non_empty_domain()`](SOMAArrayBase.html#method-non_empty_domain)
- [`tiledbsoma::SOMAArrayBase$open()`](SOMAArrayBase.html#method-open)
- [`tiledbsoma::SOMAArrayBase$print()`](SOMAArrayBase.html#method-print)
- [`tiledbsoma::SOMAArrayBase$schema()`](SOMAArrayBase.html#method-schema)
- [`tiledbsoma::SOMAArrayBase$shape()`](SOMAArrayBase.html#method-shape)
- [`tiledbsoma::SOMANDArrayBase$create()`](SOMANDArrayBase.html#method-create)
- [`tiledbsoma::SOMANDArrayBase$resize()`](SOMANDArrayBase.html#method-resize)
- [`tiledbsoma::SOMANDArrayBase$set_data_type()`](SOMANDArrayBase.html#method-set_data_type)
- [`tiledbsoma::SOMANDArrayBase$tiledbsoma_has_upgraded_shape()`](SOMANDArrayBase.html#method-tiledbsoma_has_upgraded_shape)
- [`tiledbsoma::SOMANDArrayBase$tiledbsoma_upgrade_shape()`](SOMANDArrayBase.html#method-tiledbsoma_upgrade_shape)

------------------------------------------------------------------------

### Method `read_arrow_table()`

Read as an [Arrow
table](https://arrow.apache.org/docs/r/reference/Table-class.html)
(lifecycle: maturing).

#### Usage

    SOMADenseNDArray$read_arrow_table(
      coords = NULL,
      result_order = "auto",
      log_level = "auto"
    )

#### Arguments

- `coords`:

  Optional `list` of integer vectors, one for each dimension, with a
  length equal to the number of values to read. If `NULL`, all values
  are read. List elements can be named when specifying a subset of
  dimensions.

- `result_order`:

  Optional order of read results. This can be one of either
  `"ROW_MAJOR, `"COL_MAJOR"`, or `"auto"\` (default).

- `log_level`:

  Optional logging level with default value of “`warn`”.

#### Returns

An [Arrow
table](https://arrow.apache.org/docs/r/reference/Table-class.html).

------------------------------------------------------------------------

### Method `read_dense_matrix()`

Read as a dense matrix (lifecycle: maturing).

#### Usage

    SOMADenseNDArray$read_dense_matrix(
      coords = NULL,
      result_order = "ROW_MAJOR",
      log_level = "warn"
    )

#### Arguments

- `coords`:

  Optional `list` of integer vectors, one for each dimension, with a
  length equal to the number of values to read. If `NULL`, all values
  are read. List elements can be named when specifying a subset of
  dimensions.

- `result_order`:

  Optional order of read results. This can be one of either
  `"ROW_MAJOR, `"COL_MAJOR"`, or `"auto"\` (default).

- `log_level`:

  Optional logging level with default value of “`warn`”.

#### Returns

A `matrix`.

------------------------------------------------------------------------

### Method [`write()`](https://rdrr.io/r/base/write.html)

Write matrix data to the array (lifecycle: maturing).  
  
**Note**: The `$write()` method is currently limited to writing from
two-dimensional matrices (lifecycle: maturing).

#### Usage

    SOMADenseNDArray$write(values, coords = NULL)

#### Arguments

- `values`:

  A `matrix`. Character dimension names are ignored because
  `SOMANDArray`s use integer indexing.

- `coords`:

  A `list` of integer vectors, one for each dimension, with a length
  equal to the number of values to write. If `NULL`, the default, the
  values are taken from the row and column names of `values`.

#### Returns

Invisibly returns `self`.

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    SOMADenseNDArray$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
uri <- withr::local_tempfile(pattern = "soma-dense-array")
mat <- matrix(stats::rnorm(100L ^ 2L), nrow = 100L, ncol = 100L)
mat[1:3, 1:5]
#>            [,1]       [,2]        [,3]       [,4]        [,5]
#> [1,]  1.1796642  0.5143278  0.90426912  1.3364468 -1.54103029
#> [2,] -0.2569212 -1.7513751  0.07964921 -0.8603562 -0.31031012
#> [3,] -1.0563361  0.8935975 -1.25882722  0.6665378 -0.02010818

(arr <- SOMADenseNDArrayCreate(uri, arrow::float64(), shape = dim(mat)))
#> <SOMADenseNDArray>
#>   uri: /tmp/Rtmpfr8mYm/soma-dense-array2a7e104b689f
#>   dimensions: soma_dim_0, soma_dim_1 
#>   attributes: soma_data 
arr$write(mat)
arr$close()

(arr <- SOMADenseNDArrayOpen(uri))
#> <SOMADenseNDArray>
#>   uri: /tmp/Rtmpfr8mYm/soma-dense-array2a7e104b689f
#>   dimensions: soma_dim_0, soma_dim_1 
#>   attributes: soma_data 
arr$read_arrow_table()
#> Table
#> 10000 rows x 1 columns
#> $soma_data <double not null>
```
