# The SOMA Object Base Class

Base class to implement shared functionality across the
[`SOMAArrayBase`](SOMAArrayBase.md) and
[`SOMACollectionBase`](SOMACollectionBase.md) classes. (lifecycle:
maturing)

## See also

Derived classes: [`SOMAArrayBase`](SOMAArrayBase.md),
[`SOMACollectionBase`](SOMACollectionBase.md)

## Active bindings

- `platform_config`:

  Platform configuration

- `context`:

  SOMAContext context object for TileDB operations

- `tiledbsoma_ctx`:

  SOMATileDBContext

- `tiledb_timestamp`:

  Time that object was opened at

- `uri`:

  The URI of the TileDB object

- `soma_type`:

  The SOMA object type

- `.tiledb_timestamp_range`:

  Time range for libtiledbsoma

## Methods

### Public methods

- [`SOMAObject$new()`](#method-SOMAObject-new)

- [`SOMAObject$is_open()`](#method-SOMAObject-is_open)

- [`SOMAObject$class()`](#method-SOMAObject-class)

- [`SOMAObject$mode()`](#method-SOMAObject-mode)

- [`SOMAObject$reopen()`](#method-SOMAObject-reopen)

- [`SOMAObject$exists()`](#method-SOMAObject-exists)

- [`SOMAObject$get_metadata()`](#method-SOMAObject-get_metadata)

- [`SOMAObject$set_metadata()`](#method-SOMAObject-set_metadata)

- [`SOMAObject$print()`](#method-SOMAObject-print)

- [`SOMAObject$clone()`](#method-SOMAObject-clone)

------------------------------------------------------------------------

### Method `new()`

Create a new SOMA object. (lifecycle: maturing)

#### Usage

    SOMAObject$new(
      uri,
      ...,
      platform_config = NULL,
      tiledbsoma_ctx = NULL,
      tiledb_timestamp = NULL,
      context = NULL
    )

#### Arguments

- `uri`:

  URI for the SOMA object

- `...`:

  Ignored

- `platform_config`:

  Optional platform configuration

- `tiledbsoma_ctx`:

  Optional (DEPRECATED) TileDB SOMA context

- `tiledb_timestamp`:

  Optional timestamp
  ([`POSIXct`](https://rdrr.io/r/base/DateTimeClasses.html)) to open the
  object at

- `context`:

  Optional `SOMAContext` object used for TileDB operations. If a context
  is not provided, then the default context will be used. Call
  `set_default_context` once before other SOMA operations to configure
  the default context.

------------------------------------------------------------------------

### Method `is_open()`

Determine if the object is open for reading or writing

#### Usage

    SOMAObject$is_open()

#### Returns

`TRUE` if the object is open, otherwise `FALSE`

------------------------------------------------------------------------

### Method [`class()`](https://rdrr.io/r/base/class.html)

Print the name of the R6 class

#### Usage

    SOMAObject$class()

#### Returns

The name of the R6 class

------------------------------------------------------------------------

### Method [`mode()`](https://rdrr.io/r/base/mode.html)

Get the mode of the object

#### Usage

    SOMAObject$mode()

#### Returns

The mode of the object, one of:

- “`CLOSED`”

- “`READ`”

- “`WRITE`”

- “`DELETE`”

------------------------------------------------------------------------

### Method `reopen()`

Close and reopen the TileDB object in a new mode

#### Usage

    SOMAObject$reopen(mode, tiledb_timestamp = NULL)

#### Arguments

- `mode`:

  New mode to open the object in; choose from:

  - “`READ`”

  - “`WRITE`”

  - “`DELETE`”

- `tiledb_timestamp`:

  Optional Datetime (POSIXct) with TileDB timestamp

#### Returns

Invisibly returns `self` opened in `mode`

------------------------------------------------------------------------

### Method [`exists()`](https://rdrr.io/r/base/exists.html)

Check if the object exists. (lifecycle: maturing)

#### Usage

    SOMAObject$exists()

#### Returns

`TRUE` if the object exists, otherwise `FALSE`

------------------------------------------------------------------------

### Method `get_metadata()`

Retrieve metadata. (lifecycle: maturing)

#### Usage

    SOMAObject$get_metadata(key = NULL)

#### Arguments

- `key`:

  The name of the metadata attribute to retrieve is not NULL.

#### Returns

A list of metadata values.

------------------------------------------------------------------------

### Method `set_metadata()`

Add list of metadata. (lifecycle: maturing)

#### Usage

    SOMAObject$set_metadata(metadata)

#### Arguments

- `metadata`:

  Named list of metadata to add

#### Returns

Invisibly returns `self`

------------------------------------------------------------------------

### Method [`print()`](https://rdrr.io/r/base/print.html)

Print-friendly representation of the object

#### Usage

    SOMAObject$print()

#### Returns

Invisibly returns `self`

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    SOMAObject$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
