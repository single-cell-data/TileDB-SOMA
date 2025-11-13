# SOMA Array Base Class

SOMA Array Base Class

SOMA Array Base Class

## Details

Virtual base class to add array-specific functionality to the
[`SOMAObject`](SOMAObject.md) class (lifecycle: maturing).

## See also

Derived classes: [`SOMADataFrame`](SOMADataFrame.md),
[`SOMANDArrayBase`](SOMANDArrayBase.md)

## Super class

[`tiledbsoma::SOMAObject`](SOMAObject.md) -\> `SOMAArrayBase`

## Methods

### Public methods

- [`SOMAArrayBase$open()`](#method-SOMAArrayBase-open)

- [`SOMAArrayBase$close()`](#method-SOMAArrayBase-close)

- [`SOMAArrayBase$allows_duplicates()`](#method-SOMAArrayBase-allows_duplicates)

- [`SOMAArrayBase$is_sparse()`](#method-SOMAArrayBase-is_sparse)

- [`SOMAArrayBase$schema()`](#method-SOMAArrayBase-schema)

- [`SOMAArrayBase$attributes()`](#method-SOMAArrayBase-attributes)

- [`SOMAArrayBase$attrnames()`](#method-SOMAArrayBase-attrnames)

- [`SOMAArrayBase$dimensions()`](#method-SOMAArrayBase-dimensions)

- [`SOMAArrayBase$dimnames()`](#method-SOMAArrayBase-dimnames)

- [`SOMAArrayBase$colnames()`](#method-SOMAArrayBase-colnames)

- [`SOMAArrayBase$index_column_names()`](#method-SOMAArrayBase-index_column_names)

- [`SOMAArrayBase$shape()`](#method-SOMAArrayBase-shape)

- [`SOMAArrayBase$maxshape()`](#method-SOMAArrayBase-maxshape)

- [`SOMAArrayBase$non_empty_domain()`](#method-SOMAArrayBase-non_empty_domain)

- [`SOMAArrayBase$ndim()`](#method-SOMAArrayBase-ndim)

- [`SOMAArrayBase$print()`](#method-SOMAArrayBase-print)

- [`SOMAArrayBase$clone()`](#method-SOMAArrayBase-clone)

Inherited methods

- [`tiledbsoma::SOMAObject$class()`](SOMAObject.html#method-class)
- [`tiledbsoma::SOMAObject$exists()`](SOMAObject.html#method-exists)
- [`tiledbsoma::SOMAObject$get_metadata()`](SOMAObject.html#method-get_metadata)
- [`tiledbsoma::SOMAObject$initialize()`](SOMAObject.html#method-initialize)
- [`tiledbsoma::SOMAObject$is_open()`](SOMAObject.html#method-is_open)
- [`tiledbsoma::SOMAObject$mode()`](SOMAObject.html#method-mode)
- [`tiledbsoma::SOMAObject$reopen()`](SOMAObject.html#method-reopen)
- [`tiledbsoma::SOMAObject$set_metadata()`](SOMAObject.html#method-set_metadata)

------------------------------------------------------------------------

### Method [`open()`](https://rdrr.io/r/base/connections.html)

Open the SOMA object for read or write.  
  
**Note**: [`open()`](https://rdrr.io/r/base/connections.html) is
considered internal and should not be called directly; use factory
functions (eg. [`SOMASparseNDArrayOpen()`](SOMASparseNDArrayOpen.md))
instead.

#### Usage

    SOMAArrayBase$open(mode = c("READ", "WRITE", "DELETE"))

#### Arguments

- `mode`:

  Mode to open the object in.

#### Returns

Return s`self`.

------------------------------------------------------------------------

### Method [`close()`](https://rdrr.io/r/base/connections.html)

Close the SOMA array.

#### Usage

    SOMAArrayBase$close()

#### Returns

Invisibly returns `self`.

------------------------------------------------------------------------

### Method `allows_duplicates()`

Does an array allow duplicates?

#### Usage

    SOMAArrayBase$allows_duplicates()

#### Returns

`TRUE` if the underlying TileDB array allows duplicates; otherwise
`FALSE`.

------------------------------------------------------------------------

### Method `is_sparse()`

Is an array sparse?

#### Usage

    SOMAArrayBase$is_sparse()

#### Returns

`TRUE` if the underlying TileDB array is sparse; otherwise `FALSE`.

------------------------------------------------------------------------

### Method `schema()`

Retrieve the array schema as an Arrow schema (lifecycle: maturing).

#### Usage

    SOMAArrayBase$schema()

#### Returns

An Arrow
[`Schema`](https://arrow.apache.org/docs/r/reference/Schema-class.html)
object.

------------------------------------------------------------------------

### Method [`attributes()`](https://rdrr.io/r/base/attributes.html)

Retrieve the array attributes.

#### Usage

    SOMAArrayBase$attributes()

#### Returns

A named list of array attributes; each entry contains the following
named entries:

- “`name`”: name of the attribute.

- “`type`”: datatype of the attribute.

- “`ncells`”: number of values per dimension cell.

- “`nullable`”: is the attribute nullable.

- “`filter_list`”: a list with filter information; this list contains
  the following entries:

  - “`filter_type`”

  - “`compression_level`”

  - “`bit_width`”

  - “`positive_delta`”

  - “`float_bytewidth`”

  - “`float_factor`”

  - “`float_offset`”

------------------------------------------------------------------------

### Method `attrnames()`

Retrieve attribute names (lifecycle: maturing).

#### Usage

    SOMAArrayBase$attrnames()

#### Returns

A character vector with the array's attribute names.

------------------------------------------------------------------------

### Method `dimensions()`

Retrieve the array dimensions (lifecycle: maturing)

#### Usage

    SOMAArrayBase$dimensions()

#### Returns

A named list of array dimensions; each entry contains the following
named entries:

- “`name`”: name of the dimension.

- “`type`”: datatype of the dimension.

- “`ncells`”: number of values per dimension cell.

- “`domain`”: domain of the dimension.

- “`tile`”: tile of the dimension.

- “`filter_list`”: a list with filter information; this list contains
  the following entries:

  - “`filter_type`”

  - “`compression_level`”

  - “`bit_width`”

  - “`positive_delta`”

  - “`float_bytewidth`”

  - “`float_factor`”

  - “`float_offset`”

------------------------------------------------------------------------

### Method [`dimnames()`](https://rdrr.io/r/base/dimnames.html)

Retrieve dimension names (lifecycle: maturing).

#### Usage

    SOMAArrayBase$dimnames()

#### Returns

A character vector with the array's dimension names.

------------------------------------------------------------------------

### Method [`colnames()`](https://rdrr.io/r/base/colnames.html)

Retrieve the names of all columns, including dimensions and attributes
(lifecycle: maturing).

#### Usage

    SOMAArrayBase$colnames()

#### Returns

A character vector with the array's column names.

------------------------------------------------------------------------

### Method `index_column_names()`

Retrieve names of index (dimension) columns (lifecycle: maturing)

#### Usage

    SOMAArrayBase$index_column_names()

#### Returns

A character vector with the array index (dimension) names

------------------------------------------------------------------------

### Method `shape()`

Retrieve the shape, i.e. the capacity of each dimension Attempted reads
and writes outside the `shape` will result in a run-time error: this is
the purpose of `shape`. This will not necessarily match the bounds of
occupied cells within the array. Using `$resize()`, this may be
increased up to the hard limit which `maxshape()` reports (lifecycle:
maturing).

#### Usage

    SOMAArrayBase$shape()

#### Returns

A named vector of dimension length and of the same type as the
dimension.

------------------------------------------------------------------------

### Method `maxshape()`

Retrieve the hard limit up to which the array may be resized using the
`$resize()` method (lifecycle: maturing).

#### Usage

    SOMAArrayBase$maxshape()

#### Returns

A named vector of dimension length and of the same type as the
dimension.

------------------------------------------------------------------------

### Method `non_empty_domain()`

Returns a named list of minimum/maximum pairs, one per index column,
which are the smallest and largest values written on that index column.

#### Usage

    SOMAArrayBase$non_empty_domain(index1 = FALSE, max_only = FALSE)

#### Arguments

- `index1`:

  Return the non-empty domain with 1-based indices

- `max_only`:

  Return only the max value per dimension, and return this as a vector.
  Names are dropped (lifecycle: maturing).

#### Returns

Named list of minimum/maximum values, or integer vector of maximum
values.

------------------------------------------------------------------------

### Method `ndim()`

Retrieve number of dimensions (lifecycle: maturing).

#### Usage

    SOMAArrayBase$ndim()

#### Returns

A scalar with the number of dimensions.

------------------------------------------------------------------------

### Method [`print()`](https://rdrr.io/r/base/print.html)

Print-friendly representation of the object.

#### Usage

    SOMAArrayBase$print()

#### Returns

Invisibly returns `self`.

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    SOMAArrayBase$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
