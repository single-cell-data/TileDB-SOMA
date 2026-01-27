# SOMA Sparse Nd-Array

`SOMASparseNDArray` is a sparse, N-dimensional array with offset
(zero-based) integer indexing on each dimension. The `SOMASparseNDArray`
has a user-defined schema, which includes:

- `type`: a `primitive` type, expressed as an Arrow type (e.g.,
  [`int64`](https://arrow.apache.org/docs/r/reference/data-type.html),
  [`float32`](https://arrow.apache.org/docs/r/reference/data-type.html),
  etc), indicating the type of data contained within the array.

- `shape`: the shape of the array, i.e., number and length of each
  dimension. This is a soft limit which can be increased using
  `$resize()` up to the `maxshape`.

- `maxshape`: the hard limit up to which `shape` may be increased using
  `$resize()`.

All dimensions must have a positive, non-zero length.

## Note

In TileDB this is an sparse array with `N` `int64` dimensions of domain
`[0, maxInt64)` and a single attribute.

## Duplicate Writes

As duplicate index values are not allowed, index values already present
in the object are overwritten and new index values are added (lifecycle:
maturing).

## Super classes

[`tiledbsoma::SOMAObject`](SOMAObject.md) -\>
[`tiledbsoma::SOMAArrayBase`](SOMAArrayBase.md) -\>
[`tiledbsoma::SOMANDArrayBase`](SOMANDArrayBase.md) -\>
`SOMASparseNDArray`

## Methods

### Public methods

- [`SOMASparseNDArray$read()`](#method-SOMASparseNDArray-read)

- [`SOMASparseNDArray$write()`](#method-SOMASparseNDArray-write)

- [`SOMASparseNDArray$nnz()`](#method-SOMASparseNDArray-nnz)

- [`SOMASparseNDArray$.write_coordinates()`](#method-SOMASparseNDArray-.write_coordinates)

- [`SOMASparseNDArray$clone()`](#method-SOMASparseNDArray-clone)

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

### Method `read()`

Reads a user-defined slice of the `SOMASparseNDArray`.

#### Usage

    SOMASparseNDArray$read(
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

A [SOMASparseNDArrayRead](SOMASparseNDArrayRead.md).

------------------------------------------------------------------------

### Method [`write()`](https://rdrr.io/r/base/write.html)

Write matrix-like data to the array (lifecycle: maturing).

#### Usage

    SOMASparseNDArray$write(values, bbox = NULL)

#### Arguments

- `values`:

  Any `matrix`-like object coercible to a
  [`TsparseMatrix`](https://rdrr.io/pkg/Matrix/man/TsparseMatrix-class.html).
  Character dimension names are ignored because `SOMANDArray`s use
  integer indexing.

- `bbox`:

  A vector of integers describing the upper bounds of each dimension of
  `values`. Generally should be `NULL`.

#### Returns

Invisibly returns `self`.

------------------------------------------------------------------------

### Method `nnz()`

Retrieve number of non-zero elements (lifecycle: maturing).

#### Usage

    SOMASparseNDArray$nnz()

#### Returns

A scalar with the number of non-zero elements.

------------------------------------------------------------------------

### Method `.write_coordinates()`

Write a COO table to the array.

#### Usage

    SOMASparseNDArray$.write_coordinates(values)

#### Arguments

- `values`:

  A `data.frame` or [Arrow
  table](https://arrow.apache.org/docs/r/reference/Table-class.html)
  with data in COO format; must be named with the dimension and
  attribute labels of the array.

#### Returns

Invisibly returns `self`.

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    SOMASparseNDArray$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
uri <- withr::local_tempfile(pattern = "soma-sparse-array")
mat <- Matrix::rsparsematrix(100L, 100L, 0.7, repr = "T")
mat[1:3, 1:5]
#> 3 x 5 sparse Matrix of class "dgTMatrix"
#>                                
#> [1,] 0.59 .     1.30 -1.20  .  
#> [2,] 0.30 0.55  1.20 -1.50 -2.9
#> [3,] 0.41 .    -0.37  0.59  .  

(arr <- SOMASparseNDArrayCreate(uri, arrow::float64(), shape = dim(mat)))
#> <SOMASparseNDArray>
#>   uri: /tmp/Rtmpfr8mYm/soma-sparse-array2a7e20215854
#>   dimensions: soma_dim_0, soma_dim_1 
#>   attributes: soma_data 
arr$write(mat)
arr$close()

(arr <- SOMASparseNDArrayOpen(uri))
#> <SOMASparseNDArray>
#>   uri: /tmp/Rtmpfr8mYm/soma-sparse-array2a7e20215854
#>   dimensions: soma_dim_0, soma_dim_1 
#>   attributes: soma_data 
m2 <- arr$read()$sparse_matrix()$concat()
m2[1:3, 1:5]
#> 3 x 5 sparse Matrix of class "dgTMatrix"
#>                                
#> [1,] 0.59 .     1.30 -1.20  .  
#> [2,] 0.30 0.55  1.20 -1.50 -2.9
#> [3,] 0.41 .    -0.37  0.59  .  
```
