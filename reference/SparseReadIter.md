# SOMA Read Iterator Over Sparse Matrices

`SparseReadIter` is a class that allows for iteration over a reads on
[SOMASparseNDArray](SOMASparseNDArray.md).

## Super class

[`tiledbsoma::ReadIter`](ReadIter.md) -\> `SparseReadIter`

## Methods

### Public methods

- [`SparseReadIter$new()`](#method-SparseReadIter-new)

- [`SparseReadIter$concat()`](#method-SparseReadIter-concat)

- [`SparseReadIter$clone()`](#method-SparseReadIter-clone)

Inherited methods

- [`tiledbsoma::ReadIter$read_complete()`](ReadIter.html#method-read_complete)
- [`tiledbsoma::ReadIter$read_next()`](ReadIter.html#method-read_next)

------------------------------------------------------------------------

### Method `new()`

Create (lifecycle: maturing).

#### Usage

    SparseReadIter$new(sr, shape, zero_based = FALSE)

#### Arguments

- `sr`:

  Soma reader pointer.

- `shape`:

  Shape of the full matrix.

- `zero_based`:

  Logical, if `TRUE` will make iterator for
  Matrix::[dgTMatrix-class](https://rdrr.io/pkg/Matrix/man/dgTMatrix-class.html)
  otherwise [matrixZeroBasedView](matrixZeroBasedView.md).

------------------------------------------------------------------------

### Method `concat()`

Concatenate remainder of iterator.

#### Usage

    SparseReadIter$concat()

#### Returns

[matrixZeroBasedView](matrixZeroBasedView.md) of
Matrix::[sparseMatrix](https://rdrr.io/pkg/Matrix/man/sparseMatrix.html).

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    SparseReadIter$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
dir <- withr::local_tempfile(pattern = "matrix-iter")
dir.create(dir, recursive = TRUE)
(exp <- load_dataset("soma-exp-pbmc-small", dir))
#> <SOMAExperiment>
#>   uri: /tmp/RtmpbAgXbM/matrix-iter28463c662601/soma-exp-pbmc-small
qry <- exp$axis_query("RNA")
xqry <- qry$X("data")

iter <- xqry$sparse_matrix()
stopifnot(inherits(iter, "SparseReadIter"))

while (!iter$read_complete()) {
  block <- iter$read_next()
}
```
