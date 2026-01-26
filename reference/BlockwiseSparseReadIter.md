# SOMA Blockwise Read Iterator for Sparse Matrices

Class that allows for blockwise read iteration of SOMA reads as sparse
matrices.

## Super classes

[`tiledbsoma::ReadIter`](ReadIter.md) -\>
[`tiledbsoma::BlockwiseReadIterBase`](BlockwiseReadIterBase.md) -\>
`BlockwiseSparseReadIter`

## Active bindings

- `repr`:

  Representation of the sparse matrix to return.

## Methods

### Public methods

- [`BlockwiseSparseReadIter$new()`](#method-BlockwiseSparseReadIter-new)

- [`BlockwiseSparseReadIter$concat()`](#method-BlockwiseSparseReadIter-concat)

- [`BlockwiseSparseReadIter$clone()`](#method-BlockwiseSparseReadIter-clone)

Inherited methods

- [`tiledbsoma::BlockwiseReadIterBase$read_complete()`](BlockwiseReadIterBase.html#method-read_complete)
- [`tiledbsoma::BlockwiseReadIterBase$read_next()`](BlockwiseReadIterBase.html#method-read_next)

------------------------------------------------------------------------

### Method `new()`

Create.

#### Usage

    BlockwiseSparseReadIter$new(
      sr,
      array,
      coords,
      axis,
      ...,
      repr = "T",
      reindex_disable_on_axis = NA
    )

#### Arguments

- `sr`:

  SOMA read pointer

- `array`:

  Underlying [`SOMASparseNDArray`](SOMASparseNDArray.md)

- `coords`:

  Named list of [`CoordsStrider`](CoordsStrider.md) objects; must be
  named after `array$dimnames()`

- `axis`:

  Axis to iterate over in a blockwise manner

- `...`:

  Ignored

- `repr`:

  Representation of the sparse matrix to return; choose from:

  - “`T`”: returns a
    [`TsparseMatrix`](https://rdrr.io/pkg/Matrix/man/TsparseMatrix-class.html)

  - “`R`”: returns an
    [`RsparseMatrix`](https://rdrr.io/pkg/Matrix/man/RsparseMatrix-class.html)

  - “`C`”: returns a
    [`CsparseMatrix`](https://rdrr.io/pkg/Matrix/man/CsparseMatrix-class.html)

  **Note**: passing `repr` of “`R`” or “`C`” are only available if
  re-indexing is enabled on axes `0` or `1`, respectively

- `reindex_disable_on_axis`:

  Additional axes that will not be re-indexed; the following values may
  be used as shorthands for common settings:

  - “`TRUE`”: disable re-indexing on all axes

  - “`NA`”: re-index only on `axis`, disable re-indexing on all others

  - “`FALSE`”: re-index on *all* axes, do **not** disable re-indexing

------------------------------------------------------------------------

### Method `concat()`

Concatenate the remainder of the blockwise iterator.

#### Usage

    BlockwiseSparseReadIter$concat()

#### Returns

A sparse matrix (determined by `self$repr`) with the remainder of the
iterator.

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    BlockwiseSparseReadIter$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
dir <- withr::local_tempfile(pattern = "blockwise-matrix")
dir.create(dir, recursive = TRUE)
(exp <- load_dataset("soma-exp-pbmc-small", dir))
#> <SOMAExperiment>
#>   uri: /tmp/RtmpbAgXbM/blockwise-matrix2846190372af/soma-exp-pbmc-small
qry <- exp$axis_query("RNA")
xqry <- qry$X("data")

iter <- xqry$blockwise(0L, size = 20L, reindex_disable_on_axis = TRUE)$sparse_matrix()
stopifnot(inherits(iter, "BlockwiseSparseReadIter"))

while (!iter$read_complete()) {
  block <- iter$read_next()
}
```
