# SOMA Blockwise Read Iterator Base Class

Virtual base class that allows for blockwise read iteration of SOMA
reads.

## See also

Derived classes: [`BlockwiseTableReadIter`](BlockwiseTableReadIter.md),
[`BlockwiseSparseReadIter`](BlockwiseSparseReadIter.md)

## Super class

[`tiledbsoma::ReadIter`](ReadIter.md) -\> `BlockwiseReadIterBase`

## Active bindings

- `array`:

  The underlying SOMA array.

- `axis`:

  The axis to iterate over in a blockwise fashion.

- `axes_to_reindex`:

  The axes to re-index.

- `coords`:

  A list of [`CoordsStrider`](CoordsStrider.md) objects.

- `coords_axis`:

  The [`CoordsStrider`](CoordsStrider.md) for `axis`.

- `reindex_disable_on_axis`:

  Additional axes that will not be re-indexed.

- `reindexable`:

  Shorthand to see if this iterator is poised to be re-indexed or not.

## Methods

### Public methods

- [`BlockwiseReadIterBase$new()`](#method-BlockwiseReadIterBase-new)

- [`BlockwiseReadIterBase$read_complete()`](#method-BlockwiseReadIterBase-read_complete)

- [`BlockwiseReadIterBase$read_next()`](#method-BlockwiseReadIterBase-read_next)

- [`BlockwiseReadIterBase$clone()`](#method-BlockwiseReadIterBase-clone)

Inherited methods

- [`tiledbsoma::ReadIter$concat()`](ReadIter.html#method-concat)

------------------------------------------------------------------------

### Method `new()`

Create.

#### Usage

    BlockwiseReadIterBase$new(
      sr,
      array,
      coords,
      axis,
      ...,
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

- `reindex_disable_on_axis`:

  Additional axes that will not be re-indexed; the following values may
  be used as shorthands for common settings:

  - “`TRUE`”: disable re-indexing on all axes

  - “`NA`”: re-index only on `axis`, disable re-indexing on all others

  - “`FALSE`”: re-index on *all* axes, do **not** disable re-indexing

------------------------------------------------------------------------

### Method `read_complete()`

Check if the iterated read is complete or not.

#### Usage

    BlockwiseReadIterBase$read_complete()

#### Returns

`TRUE` if read is complete, otherwise `FALSE`.

------------------------------------------------------------------------

### Method `read_next()`

Read the next chunk of the iterated read. If read is complete, throws an
`iterationCompleteWarning` warning and returns `NULL`.

#### Usage

    BlockwiseReadIterBase$read_next()

#### Returns

`NULL` or the next blockwise chunk of the iterated read.

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    BlockwiseReadIterBase$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
