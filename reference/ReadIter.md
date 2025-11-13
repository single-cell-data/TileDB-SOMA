# SOMA Read Iterator Base class

SOMA Read Iterator Base class

SOMA Read Iterator Base class

## Details

Virtual class that allows for read iteration of SOMA reads.

## See also

[`BlockwiseReadIterBase`](BlockwiseReadIterBase.md),
[`SparseReadIter`](SparseReadIter.md),
[`TableReadIter`](TableReadIter.md)

## Methods

### Public methods

- [`ReadIter$new()`](#method-ReadIter-new)

- [`ReadIter$read_complete()`](#method-ReadIter-read_complete)

- [`ReadIter$read_next()`](#method-ReadIter-read_next)

- [`ReadIter$concat()`](#method-ReadIter-concat)

- [`ReadIter$clone()`](#method-ReadIter-clone)

------------------------------------------------------------------------

### Method `new()`

Create (lifecycle: maturing).

#### Usage

    ReadIter$new(sr)

#### Arguments

- `sr`:

  soma read pointer.

------------------------------------------------------------------------

### Method `read_complete()`

@description Check if iterated read is complete or not (lifecycle:
maturing).

#### Usage

    ReadIter$read_complete()

#### Returns

logical

------------------------------------------------------------------------

### Method `read_next()`

Read the next chunk of an iterated read. If read is complete, returns
`NULL` and raises warning (lifecycle: maturing).

#### Usage

    ReadIter$read_next()

#### Returns

`NULL` or one of
[arrow::Table](https://arrow.apache.org/docs/r/reference/Table-class.html),
[matrixZeroBasedView](matrixZeroBasedView.md).

------------------------------------------------------------------------

### Method `concat()`

Concatenate remainder of iterator.

#### Usage

    ReadIter$concat()

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    ReadIter$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
