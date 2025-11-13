# Blockwise Sparse ND-Array Reader

Blockwise reader for [`SOMASparseNDArray`](SOMASparseNDArray.md)

## Super class

[`tiledbsoma::SOMASparseNDArrayReadBase`](SOMASparseNDArrayReadBase.md)
-\> `SOMASparseNDArrayBlockwiseRead`

## Active bindings

- `axis`:

  The axis to iterate over in a blockwise fashion

## Methods

### Public methods

- [`SOMASparseNDArrayBlockwiseRead$new()`](#method-SOMASparseNDArrayBlockwiseRead-new)

- [`SOMASparseNDArrayBlockwiseRead$tables()`](#method-SOMASparseNDArrayBlockwiseRead-tables)

- [`SOMASparseNDArrayBlockwiseRead$sparse_matrix()`](#method-SOMASparseNDArrayBlockwiseRead-sparse_matrix)

------------------------------------------------------------------------

### Method `new()`

Create

#### Usage

    SOMASparseNDArrayBlockwiseRead$new(
      sr,
      array,
      coords,
      axis,
      ...,
      size,
      reindex_disable_on_axis = NA
    )

#### Arguments

- `sr`:

  SOMA read pointer

- `array`:

  Underlying [`SOMASparseNDArray`](SOMASparseNDArray.md)

- `coords`:

  Optional named list of
  [`integer64`](https://rdrr.io/pkg/bit64/man/bit64-package.html)
  values; must be named after `array$dimnames()`

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

------------------------------------------------------------------------

### Method `tables()`

Read as an
[`Arrow::Table`](https://arrow.apache.org/docs/r/reference/Table-class.html)

#### Usage

    SOMASparseNDArrayBlockwiseRead$tables()

#### Returns

A blockwise iterator yielding chunks as
[`Arrow::Table`](https://arrow.apache.org/docs/r/reference/Table-class.html)`s`

------------------------------------------------------------------------

### Method `sparse_matrix()`

Read as a sparse matrix

#### Usage

    SOMASparseNDArrayBlockwiseRead$sparse_matrix(repr = "T")

#### Arguments

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

#### Returns

A blockwise iterator yielding chunks as sparse matrices
