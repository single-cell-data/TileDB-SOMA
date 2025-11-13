# Zero-based Wrapper for Sparse Matrices

Zero-based Wrapper for Sparse Matrices

Zero-based Wrapper for Sparse Matrices

## Details

`matrixZeroBasedView` is a wrapper shim for a matrix or
[`Matrix::sparseMatrix`](https://rdrr.io/pkg/Matrix/man/sparseMatrix.html)
that allows elemental matrix access using zero-based indices.

## Methods

### Public methods

- [`matrixZeroBasedView$new()`](#method-matrixZeroBasedView-new)

- [`matrixZeroBasedView$take()`](#method-matrixZeroBasedView-take)

- [`matrixZeroBasedView$dim()`](#method-matrixZeroBasedView-dim)

- [`matrixZeroBasedView$nrow()`](#method-matrixZeroBasedView-nrow)

- [`matrixZeroBasedView$ncol()`](#method-matrixZeroBasedView-ncol)

- [`matrixZeroBasedView$get_one_based_matrix()`](#method-matrixZeroBasedView-get_one_based_matrix)

- [`matrixZeroBasedView$sum()`](#method-matrixZeroBasedView-sum)

- [`matrixZeroBasedView$print()`](#method-matrixZeroBasedView-print)

- [`matrixZeroBasedView$clone()`](#method-matrixZeroBasedView-clone)

------------------------------------------------------------------------

### Method `new()`

Initialize (lifecycle: maturing).

#### Usage

    matrixZeroBasedView$new(x)

#### Arguments

- `x`:

  A matrix.

------------------------------------------------------------------------

### Method `take()`

Zero-based matrix element access.

#### Usage

    matrixZeroBasedView$take(i = NULL, j = NULL)

#### Arguments

- `i`:

  Row index (zero-based).

- `j`:

  Column index (zero-based).

#### Returns

The specified matrix slice as another `matrixZeroBasedView`.

------------------------------------------------------------------------

### Method [`dim()`](https://rdrr.io/r/base/dim.html)

dim.

#### Usage

    matrixZeroBasedView$dim()

#### Returns

The dimensions of the matrix.

------------------------------------------------------------------------

### Method [`nrow()`](https://rdrr.io/r/base/nrow.html)

nrow.

#### Usage

    matrixZeroBasedView$nrow()

#### Returns

Matrix row count.

------------------------------------------------------------------------

### Method [`ncol()`](https://rdrr.io/r/base/nrow.html)

ncol.

#### Usage

    matrixZeroBasedView$ncol()

#### Returns

Matrix column count.

------------------------------------------------------------------------

### Method `get_one_based_matrix()`

Get the one-based R matrix with its original class.

#### Usage

    matrixZeroBasedView$get_one_based_matrix()

#### Returns

One-based matrix.

------------------------------------------------------------------------

### Method [`sum()`](https://rdrr.io/r/base/sum.html)

Perform arithmetic sum between this `matrixZeroBasedView` and another
`matrixZeroBasedView`.

#### Usage

    matrixZeroBasedView$sum(x)

#### Arguments

- `x`:

  the `matrixZeroBasedView` to sum.

#### Returns

The result of the sum as a `matrixZeroBasedView`.

------------------------------------------------------------------------

### Method [`print()`](https://rdrr.io/r/base/print.html)

print.

#### Usage

    matrixZeroBasedView$print()

#### Returns

Invisibly returns `self`.

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    matrixZeroBasedView$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
(mat <- Matrix::rsparsematrix(3L, 3L, 0.3))
#> 3 x 3 sparse Matrix of class "dgCMatrix"
#>                  
#> [1,]  .    0.84 .
#> [2,] -0.27 .    .
#> [3,]  .    0.24 .
(mat0 <- matrixZeroBasedView$new(mat))
#> Non-mutable 0-based 'view' class for matrices.
#> To get 1-based matrix use `x$get_one_based_matrix()
#> Dimensions: 3x3

mat0$take(0, 0)
#> Non-mutable 0-based 'view' class for matrices.
#> To get 1-based matrix use `x$get_one_based_matrix()
#> Dimensions: 1x1
mat0$take(0, 0:2)$get_one_based_matrix()
#> 1 x 3 sparse Matrix of class "dgCMatrix"
#>              
#> [1,] . 0.84 .
```
