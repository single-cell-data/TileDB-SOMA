# SOMA ND-Array Base Class

Virtual base class to add ND-array-specific functionality to the
[`SOMAArrayBase`](SOMAArrayBase.md) class (lifecycle: maturing).

## Lifecycle

As of tiledbsoma 2.1.0, `$set_data_type()` is deprecated; this
functionality is no longer needed as `libtiledbsoma` now accurately sets
the TileDB data type upon array opening

## See also

Derived classes: [`SOMADenseNDArray`](SOMADenseNDArray.md),
[`SOMASparseNDArray`](SOMASparseNDArray.md)

## Super classes

[`tiledbsoma::SOMAObject`](SOMAObject.md) -\>
[`tiledbsoma::SOMAArrayBase`](SOMAArrayBase.md) -\> `SOMANDArrayBase`

## Methods

### Public methods

- [`SOMANDArrayBase$create()`](#method-SOMANDArrayBase-create)

- [`SOMANDArrayBase$set_data_type()`](#method-SOMANDArrayBase-set_data_type)

- [`SOMANDArrayBase$tiledbsoma_has_upgraded_shape()`](#method-SOMANDArrayBase-tiledbsoma_has_upgraded_shape)

- [`SOMANDArrayBase$resize()`](#method-SOMANDArrayBase-resize)

- [`SOMANDArrayBase$tiledbsoma_upgrade_shape()`](#method-SOMANDArrayBase-tiledbsoma_upgrade_shape)

- [`SOMANDArrayBase$clone()`](#method-SOMANDArrayBase-clone)

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

------------------------------------------------------------------------

### Method `create()`

Create a SOMA NDArray named with the URI (lifecycle: maturing).  
  
**Note**: `$create()` is considered internal and should not be called
directly; use factory functions (eg.
[`SOMASparseNDArrayCreate()`](SOMASparseNDArrayCreate.md)) instead.

#### Usage

    SOMANDArrayBase$create(type, shape, platform_config = NULL)

#### Arguments

- `type`:

  An [Arrow
  type](https://arrow.apache.org/docs/r/reference/data-type.html)
  defining the type of each element in the array.

- `shape`:

  a vector of integers defining the shape of the array.

- `platform_config`:

  A [platform configuration](PlatformConfig.md) object

#### Returns

Returns `self`.

------------------------------------------------------------------------

### Method `set_data_type()`

Sets a cache value for the datatype.

#### Usage

    SOMANDArrayBase$set_data_type(type)

#### Arguments

- `type`:

  A character value describing the TileDB data type.

------------------------------------------------------------------------

### Method `tiledbsoma_has_upgraded_shape()`

Test if the array has the upgraded resizeable shape feature from
TileDB-SOMA 1.15, the array was created with this support, or it has had
`$upgrade_domain()` applied to it (lifecycle: maturing).

#### Usage

    SOMANDArrayBase$tiledbsoma_has_upgraded_shape()

#### Returns

Returns `TRUE` if the array has the upgraded resizable shape feature;
otherwise, returns `FALSE`.

------------------------------------------------------------------------

### Method `resize()`

Increases the shape of the array as specified, up to the hard limit
which is `maxshape`. Raises an error if the new shape is less than the
current shape or exceeds `maxshape` in any dimension. Also raises an
error if the array doesn't already have a shape; in that case please
call `$tiledbsoma_upgrade_shape()` (lifecycle: maturing).

#### Usage

    SOMANDArrayBase$resize(new_shape, check_only = FALSE)

#### Arguments

- `new_shape`:

  An integerish vector of the same length as the array's `$ndim()`.

- `check_only`:

  If true, does not apply the operation, but only reports whether it
  would have succeeded.

#### Returns

If `check_only`, returns the empty string if no error is detected, else
a description of the error. Otherwise, invisibly returns `NULL`.

------------------------------------------------------------------------

### Method `tiledbsoma_upgrade_shape()`

Allows the array to have a resizeable shape as described in the
TileDB-SOMA 1.15 release notes. Raises an error if the shape exceeds
`maxshape` in any dimension, or if the array already has a shape. The
methods `$tiledbsoma_upgrade_shape()` and `$resize()` are very similar:
the former must be called on a pre-1.15 array the first time a shape is
set on it; the latter must be used for subsequent resizes on any array
which already has upgraded shape (lifecycle: maturing).

#### Usage

    SOMANDArrayBase$tiledbsoma_upgrade_shape(shape, check_only = FALSE)

#### Arguments

- `shape`:

  An integerish vector of the same length as the array's `$ndim()`.

- `check_only`:

  If true, does not apply the operation, but only reports whether it
  would have succeeded.

#### Returns

If `check_only`, returns the empty string if no error is detected, else
a description of the error. Otherwise, invisibly returns `NULL`.

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    SOMANDArrayBase$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
