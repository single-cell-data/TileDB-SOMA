# SOMA Sparse ND-Array Reader Base

Base class for SOMA sparse ND-array reads

## Active bindings

- `sr`:

  The SOMA read pointer

- `array`:

  The underlying [`SOMASparseNDArray`](SOMASparseNDArray.md)

- `coords`:

  The iterated coordinates for the read

- `shape`:

  The shape of the underlying array

## Methods

### Public methods

- [`SOMASparseNDArrayReadBase$new()`](#method-SOMASparseNDArrayReadBase-new)

------------------------------------------------------------------------

### Method `new()`

Create

#### Usage

    SOMASparseNDArrayReadBase$new(sr, array, coords = NULL)

#### Arguments

- `sr`:

  SOMA read pointer

- `array`:

  Underlying [`SOMASparseNDArray`](SOMASparseNDArray.md)

- `coords`:

  Optional named list of
  [`integer64`](https://rdrr.io/pkg/bit64/man/bit64-package.html)
  values; must be named after `array$dimnames()`
