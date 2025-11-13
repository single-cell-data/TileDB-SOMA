# Ephemeral Collection Base

Virtual base class for ephemeral collections; ephemeral collections are
equivalent to [SOMA collections](SOMACollection.md) but are stored in
memory instead of on disk.

## See also

Derived classes: [`EphemeralCollection`](EphemeralCollection.md),
[`EphemeralMeasurement`](EphemeralMeasurement.md),
[`EphemeralExperiment`](EphemeralExperiment.md)

## Super classes

[`tiledbsoma::SOMAObject`](SOMAObject.md) -\>
[`tiledbsoma::SOMACollectionBase`](SOMACollectionBase.md) -\>
`EphemeralCollectionBase`

## Active bindings

- `uri`:

  “`ephemeral-collection:<MEMORY_ADDRESS>`”.

- `members`:

  A list with the members of this collection

- `soma_type`:

  Dummy field for ephemeral objects for compatibility with SOMA
  collections.

- `platform_config`:

  Dummy field for ephemeral objects for compatibility with SOMA
  collections.

- `tiledbsoma_ctx`:

  Dummy field for ephemeral objects for compatibility with SOMA
  collections.

- `tiledb_timestamp`:

  Dummy field for ephemeral objects for compatibility with SOMA
  collections

- `.tiledb_timestamp_range`:

  Dummy field for ephemeral objects for compatibility with SOMA
  collections

## Methods

### Public methods

- [`EphemeralCollectionBase$new()`](#method-EphemeralCollectionBase-new)

- [`EphemeralCollectionBase$create()`](#method-EphemeralCollectionBase-create)

- [`EphemeralCollectionBase$open()`](#method-EphemeralCollectionBase-open)

- [`EphemeralCollectionBase$close()`](#method-EphemeralCollectionBase-close)

- [`EphemeralCollectionBase$exists()`](#method-EphemeralCollectionBase-exists)

- [`EphemeralCollectionBase$print()`](#method-EphemeralCollectionBase-print)

- [`EphemeralCollectionBase$get_tiledb_config()`](#method-EphemeralCollectionBase-get_tiledb_config)

- [`EphemeralCollectionBase$length()`](#method-EphemeralCollectionBase-length)

- [`EphemeralCollectionBase$names()`](#method-EphemeralCollectionBase-names)

- [`EphemeralCollectionBase$set()`](#method-EphemeralCollectionBase-set)

- [`EphemeralCollectionBase$get()`](#method-EphemeralCollectionBase-get)

- [`EphemeralCollectionBase$remove()`](#method-EphemeralCollectionBase-remove)

- [`EphemeralCollectionBase$set_metadata()`](#method-EphemeralCollectionBase-set_metadata)

- [`EphemeralCollectionBase$get_metadata()`](#method-EphemeralCollectionBase-get_metadata)

- [`EphemeralCollectionBase$add_new_collection()`](#method-EphemeralCollectionBase-add_new_collection)

- [`EphemeralCollectionBase$add_new_dataframe()`](#method-EphemeralCollectionBase-add_new_dataframe)

- [`EphemeralCollectionBase$add_new_dense_ndarray()`](#method-EphemeralCollectionBase-add_new_dense_ndarray)

- [`EphemeralCollectionBase$add_new_sparse_ndarray()`](#method-EphemeralCollectionBase-add_new_sparse_ndarray)

- [`EphemeralCollectionBase$clone()`](#method-EphemeralCollectionBase-clone)

Inherited methods

- [`tiledbsoma::SOMAObject$class()`](SOMAObject.html#method-class)
- [`tiledbsoma::SOMAObject$is_open()`](SOMAObject.html#method-is_open)
- [`tiledbsoma::SOMAObject$mode()`](SOMAObject.html#method-mode)
- [`tiledbsoma::SOMAObject$reopen()`](SOMAObject.html#method-reopen)

------------------------------------------------------------------------

### Method `new()`

Create an ephemeral collection.

#### Usage

    EphemeralCollectionBase$new(...)

#### Arguments

- `...`:

  Ignored

------------------------------------------------------------------------

### Method `create()`

Create a new, empty ephemeral collection.

#### Usage

    EphemeralCollectionBase$create()

#### Returns

Returns a new ephemeral collection of class `class(self)`.

------------------------------------------------------------------------

### Method [`open()`](https://rdrr.io/r/base/connections.html)

Dummy method for ephemeral objects for compatibility with SOMA
collections.

#### Usage

    EphemeralCollectionBase$open(mode)

#### Arguments

- `mode`:

  Ignored for ephemeral objects.

#### Returns

Throws an error as this method is not supported by ephemeral objects.

------------------------------------------------------------------------

### Method [`close()`](https://rdrr.io/r/base/connections.html)

Dummy method for ephemeral objects for compatibility with SOMA
collections.

#### Usage

    EphemeralCollectionBase$close()

#### Returns

Invisibly returns `NULL`.

------------------------------------------------------------------------

### Method [`exists()`](https://rdrr.io/r/base/exists.html)

Dummy method for ephemeral objects for compatibility with SOMA
collections.

#### Usage

    EphemeralCollectionBase$exists()

#### Returns

Returns `FALSE` as ephemeral collections do not exist on-disk.

------------------------------------------------------------------------

### Method [`print()`](https://rdrr.io/r/base/print.html)

Special method for printing object representation to console.

#### Usage

    EphemeralCollectionBase$print()

#### Returns

Prints details about the ephemeral collection and invisibly returns
itself.

------------------------------------------------------------------------

### Method `get_tiledb_config()`

Dummy method for ephemeral objects for compatibility with SOMA
collections.

#### Usage

    EphemeralCollectionBase$get_tiledb_config(param = NULL)

#### Arguments

- `param`:

  Ignored for ephemeral objects.

#### Returns

Returns `NULL` as ephemeral collections do not have an on-disk
configuration.

------------------------------------------------------------------------

### Method [`length()`](https://rdrr.io/r/base/length.html)

Retrieve the number of items in the collection.

#### Usage

    EphemeralCollectionBase$length()

#### Returns

The length of the collection.

------------------------------------------------------------------------

### Method [`names()`](https://rdrr.io/r/base/names.html)

Retrieve the names of members. (lifecycle: maturing).

#### Usage

    EphemeralCollectionBase$names()

#### Returns

A `character` vector of member names.

------------------------------------------------------------------------

### Method `set()`

Add object to an ephemeral collection.

#### Usage

    EphemeralCollectionBase$set(object, name = NULL, relative = NULL)

#### Arguments

- `object`:

  A SOMA object (eg. [`SOMACollection`](SOMACollection.md)) to add to
  the collection.

- `name`:

  A name to add `object` as.

- `relative`:

  Ignored for ephemeral objects.

#### Returns

\[chainable\] Invisibly returns `self` with `object` added as `name`.

------------------------------------------------------------------------

### Method [`get()`](https://rdrr.io/r/base/get.html)

Get objects from an ephemeral collection.

#### Usage

    EphemeralCollectionBase$get(name)

#### Arguments

- `name`:

  Name of object in the collection to get.

#### Returns

The object named `name`.

------------------------------------------------------------------------

### Method [`remove()`](https://rdrr.io/r/base/rm.html)

Remove objects from an ephemeral collection.

#### Usage

    EphemeralCollectionBase$remove(name)

#### Arguments

- `name`:

  Name of object to remove from the collection.

#### Returns

\[chainable\] Invisibly returns `self` with the object at `name`
removed.

------------------------------------------------------------------------

### Method `set_metadata()`

Dummy method for ephemeral objects for compatibility with SOMA
collections.

#### Usage

    EphemeralCollectionBase$set_metadata(metadata)

#### Arguments

- `metadata`:

  Ignored for ephemeral objects.

#### Returns

Throws an error as this method is not supported by ephemeral objects.

------------------------------------------------------------------------

### Method `get_metadata()`

Dummy method for ephemeral objects for compatibility with SOMA
collections.

#### Usage

    EphemeralCollectionBase$get_metadata(key = NULL)

#### Arguments

- `key`:

  Ignored for ephemeral objects.

#### Returns

An empty list.

------------------------------------------------------------------------

### Method `add_new_collection()`

Dummy method for ephemeral objects for compatibility with SOMA
collections.

#### Usage

    EphemeralCollectionBase$add_new_collection(object, key)

#### Arguments

- `object, key`:

  Ignored for ephemeral objects.

#### Returns

Throws an error as this method is not supported by ephemeral objects.

------------------------------------------------------------------------

### Method `add_new_dataframe()`

Dummy method for ephemeral objects for compatibility with SOMA
collections.

#### Usage

    EphemeralCollectionBase$add_new_dataframe(key, schema, index_column_names)

#### Arguments

- `key, schema, index_column_names`:

  Ignored for ephemeral objects.

#### Returns

Throws an error as this method is not supported by ephemeral objects.

------------------------------------------------------------------------

### Method `add_new_dense_ndarray()`

Dummy method for ephemeral objects for compatibility with SOMA
collections.

#### Usage

    EphemeralCollectionBase$add_new_dense_ndarray(key, type, shape)

#### Arguments

- `key, type, shape`:

  Ignored for ephemeral objects.

#### Returns

Throws an error as this method is not supported by ephemeral objects.

------------------------------------------------------------------------

### Method `add_new_sparse_ndarray()`

Dummy method for ephemeral objects for compatibility with SOMA
collections.

#### Usage

    EphemeralCollectionBase$add_new_sparse_ndarray(key, type, shape)

#### Arguments

- `key, type, shape`:

  Ignored for ephemeral objects.

#### Returns

Throws an error as this method is not supported by ephemeral objects.

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    EphemeralCollectionBase$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
