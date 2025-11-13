# SOMASparseNDArrayRead

Intermediate type to choose result format when reading a sparse array

## Super class

[`tiledbsoma::SOMASparseNDArrayReadBase`](SOMASparseNDArrayReadBase.md)
-\> `SOMASparseNDArrayRead`

## Methods

### Public methods

- [`SOMASparseNDArrayRead$sparse_matrix()`](#method-SOMASparseNDArrayRead-sparse_matrix)

- [`SOMASparseNDArrayRead$tables()`](#method-SOMASparseNDArrayRead-tables)

- [`SOMASparseNDArrayRead$blockwise()`](#method-SOMASparseNDArrayRead-blockwise)

Inherited methods

- [`tiledbsoma::SOMASparseNDArrayReadBase$initialize()`](SOMASparseNDArrayReadBase.html#method-initialize)

------------------------------------------------------------------------

### Method `sparse_matrix()`

Read as a sparse matrix (lifecycle: maturing). Returns an iterator of
Matrix::[dgTMatrix-class](https://rdrr.io/pkg/Matrix/man/dgTMatrix-class.html)
or [matrixZeroBasedView](matrixZeroBasedView.md) of it.

#### Usage

    SOMASparseNDArrayRead$sparse_matrix(zero_based = FALSE)

#### Arguments

- `zero_based`:

  Logical, if `TRUE` returns iterator of
  [matrixZeroBasedView](matrixZeroBasedView.md) if `FALSE` returns
  iterator of
  Matrix::[dgTMatrix-class](https://rdrr.io/pkg/Matrix/man/dgTMatrix-class.html).

#### Returns

[SparseReadIter](SparseReadIter.md)

------------------------------------------------------------------------

### Method `tables()`

Read as a
arrow::[Table](https://arrow.apache.org/docs/r/reference/Table-class.html)
(lifecycle: maturing). Returns an iterator of
arrow::[Table](https://arrow.apache.org/docs/r/reference/Table-class.html).

#### Usage

    SOMASparseNDArrayRead$tables()

#### Returns

[TableReadIter](TableReadIter.md)

------------------------------------------------------------------------

### Method `blockwise()`

Read in a blockwise fashion

#### Usage

    SOMASparseNDArrayRead$blockwise(
      axis,
      ...,
      size = NULL,
      reindex_disable_on_axis = NA
    )

#### Arguments

- `axis`:

  Axis to iterate over in a blockwise manner

- `...`:

  Ignored

- `size`:

  The size of each blockwise chunk to generate

- `reindex_disable_on_axis`:

  Additional axes that will not be re-indexed; the following values may
  be used as shorthands for common settings:

  - “`TRUE`”: disable re-indexing on all axes

  - “`NA`”: re-index only on `axis`, disable re-indexing on all others

  - “`FALSE`”: re-index on *all* axes, do **not** disable re-indexing

#### Returns

A [`SOMASparseNDArrayBlockwiseRead`](SOMASparseNDArrayBlockwiseRead.md)
iterated reader
