# SOMA Collection Base Class

Base class for objects containing persistent collection of SOMA objects,
mapping string keys to any SOMA object (lifecycle: maturing).

## Note

This method is not supported for Carrara (TileDB v3) URIs. For Carrara
collections, use the `add_new_*` methods instead, which create child
objects at nested URIs that are automatically registered with the parent
collection.

## Carrara (TileDB v3) behavior

For Carrara URIs, child objects created at nested URIs are automatically
added to the parent collection. Calling this method on an already-
registered child is a **no-op** for backward compatibility.

## See also

Derived classes: [`SOMACollection`](SOMACollection.md),
[`SOMAMeasurement`](SOMAMeasurement.md),
[`SOMAExperiment`](SOMAExperiment.md)

## Super class

[`tiledbsoma::SOMAObject`](SOMAObject.md) -\> `SOMACollectionBase`

## Active bindings

- `members`:

  A list with the members of this collection.

## Methods

### Public methods

- [`SOMACollectionBase$create()`](#method-SOMACollectionBase-create)

- [`SOMACollectionBase$open()`](#method-SOMACollectionBase-open)

- [`SOMACollectionBase$close()`](#method-SOMACollectionBase-close)

- [`SOMACollectionBase$set()`](#method-SOMACollectionBase-set)

- [`SOMACollectionBase$get()`](#method-SOMACollectionBase-get)

- [`SOMACollectionBase$remove()`](#method-SOMACollectionBase-remove)

- [`SOMACollectionBase$length()`](#method-SOMACollectionBase-length)

- [`SOMACollectionBase$names()`](#method-SOMACollectionBase-names)

- [`SOMACollectionBase$set_metadata()`](#method-SOMACollectionBase-set_metadata)

- [`SOMACollectionBase$add_new_collection()`](#method-SOMACollectionBase-add_new_collection)

- [`SOMACollectionBase$add_new_dataframe()`](#method-SOMACollectionBase-add_new_dataframe)

- [`SOMACollectionBase$add_new_dense_ndarray()`](#method-SOMACollectionBase-add_new_dense_ndarray)

- [`SOMACollectionBase$add_new_sparse_ndarray()`](#method-SOMACollectionBase-add_new_sparse_ndarray)

- [`SOMACollectionBase$print()`](#method-SOMACollectionBase-print)

- [`SOMACollectionBase$clone()`](#method-SOMACollectionBase-clone)

Inherited methods

- [`tiledbsoma::SOMAObject$class()`](SOMAObject.html#method-class)
- [`tiledbsoma::SOMAObject$exists()`](SOMAObject.html#method-exists)
- [`tiledbsoma::SOMAObject$get_metadata()`](SOMAObject.html#method-get_metadata)
- [`tiledbsoma::SOMAObject$initialize()`](SOMAObject.html#method-initialize)
- [`tiledbsoma::SOMAObject$is_open()`](SOMAObject.html#method-is_open)
- [`tiledbsoma::SOMAObject$mode()`](SOMAObject.html#method-mode)
- [`tiledbsoma::SOMAObject$reopen()`](SOMAObject.html#method-reopen)

------------------------------------------------------------------------

### Method `create()`

Create a SOMA collection (lifecycle: maturing).  
  
**Note**: `$create()` is considered internal and should not be called
directly; use factory functions (eg.
[`SOMACollectionCreate()`](SOMACollectionCreate.md)) instead.

#### Usage

    SOMACollectionBase$create()

#### Returns

Returns `self`.

------------------------------------------------------------------------

### Method [`open()`](https://rdrr.io/r/base/connections.html)

Open the SOMA collection for read or write.  
  
**Note**: `$open()` is considered internal and should not be called#'
directly; use factory functions (eg.
[`SOMACollectionOpen()`](SOMACollectionOpen.md)) instead.

#### Usage

    SOMACollectionBase$open(mode = c("READ", "WRITE", "DELETE"))

#### Arguments

- `mode`:

  Mode to open the object in.

#### Returns

Returns `self`.

------------------------------------------------------------------------

### Method [`close()`](https://rdrr.io/r/base/connections.html)

Close the SOMA collection.

#### Usage

    SOMACollectionBase$close()

#### Returns

Invisibly returns `self`.

------------------------------------------------------------------------

### Method `set()`

Add a new SOMA object to the collection (lifecycle: maturing).

This method adds an existing SOMA object to the collection under the
specified key. Replacing an existing key is not supported; attempting to
add an object with a key that already exists will raise an error.

#### Usage

    SOMACollectionBase$set(object, name = NULL, relative = NULL)

#### Arguments

- `object`:

  SOMA object.

- `name`:

  The name to use for the object; defaults to the basename of
  `object$uri`.

- `relative`:

  An optional logical value indicating whether the new object's URI is
  relative to the collection's URI. If `NULL` (the default), the
  object's URI is assumed to be relative unless it is a `tiledb://` URI.

#### Returns

Invisibly returns `self`.

------------------------------------------------------------------------

### Method [`get()`](https://rdrr.io/r/base/get.html)

Retrieve a SOMA object by name. (lifecycle: maturing)

#### Usage

    SOMACollectionBase$get(name)

#### Arguments

- `name`:

  The name of the object to retrieve.

#### Returns

The SOMA object stored as `name`.

------------------------------------------------------------------------

### Method [`remove()`](https://rdrr.io/r/base/rm.html)

Remove member. (lifecycle: maturing)

#### Usage

    SOMACollectionBase$remove(name)

#### Arguments

- `name`:

  Name of the member to remove.

#### Returns

Invisibly returns `self`

------------------------------------------------------------------------

### Method [`length()`](https://rdrr.io/r/base/length.html)

Get the number of members in this collection (lifecycle: maturing).

#### Usage

    SOMACollectionBase$length()

#### Returns

The number of members in this collection.

------------------------------------------------------------------------

### Method [`names()`](https://rdrr.io/r/base/names.html)

Retrieve the names of members (lifecycle: maturing).

#### Usage

    SOMACollectionBase$names()

#### Returns

A character vector of member names.

------------------------------------------------------------------------

### Method `set_metadata()`

Add list of metadata (lifecycle: maturing).

#### Usage

    SOMACollectionBase$set_metadata(metadata)

#### Arguments

- `metadata`:

  Named list of metadata to add.

#### Returns

Invisibly returns `self`.

------------------------------------------------------------------------

### Method `add_new_collection()`

Add a new SOMA collection to this collection (lifecycle: maturing).

#### Usage

    SOMACollectionBase$add_new_collection(object, key)

#### Arguments

- `object`:

  SOMA collection object.

- `key`:

  The key to be added.

#### Returns

Returns `object`.

------------------------------------------------------------------------

### Method `add_new_dataframe()`

Add a new SOMA data frame to this collection (lifecycle: maturing).

#### Usage

    SOMACollectionBase$add_new_dataframe(
      key,
      schema,
      index_column_names,
      domain,
      platform_config = NULL
    )

#### Arguments

- `key`:

  The key to be added.

- `schema`:

  Arrow schema argument passed on to `SOMADataFrame$create()`.

- `index_column_names`:

  Index column names passed on to `SOMADataFrame$create()`.

- `domain`:

  As in [`SOMADataFrameCreate`](SOMADataFrameCreate.md).

- `platform_config`:

  A [platform configuration](PlatformConfig.md) object

- `platform_config`:

  A [platform configuration](PlatformConfig.md) object

- `platform_config`:

  A [platform configuration](PlatformConfig.md) object

#### Returns

Returns the newly created data frame stored at `key`.

------------------------------------------------------------------------

### Method `add_new_dense_ndarray()`

Add a new SOMA DenseNdArray to this collection (lifecycle: maturing).

#### Usage

    SOMACollectionBase$add_new_dense_ndarray(
      key,
      type,
      shape,
      platform_config = NULL
    )

#### Arguments

- `key`:

  The key to be added.

- `type`:

  An [Arrow
  type](https://arrow.apache.org/docs/r/reference/data-type.html)
  defining the type of each element in the array.

- `shape`:

  a vector of integers defining the shape of the array.

- `platform_config`:

  A [platform configuration](PlatformConfig.md) object

- `platform_config`:

  A [platform configuration](PlatformConfig.md) object

- `platform_config`:

  A [platform configuration](PlatformConfig.md) object

#### Returns

Returns the newly-created array stored at `key`.

------------------------------------------------------------------------

### Method `add_new_sparse_ndarray()`

Add a new SOMA SparseNdArray to this collection (lifecycle: maturing).

#### Usage

    SOMACollectionBase$add_new_sparse_ndarray(
      key,
      type,
      shape,
      platform_config = NULL
    )

#### Arguments

- `key`:

  The key to be added.

- `type`:

  An [Arrow
  type](https://arrow.apache.org/docs/r/reference/data-type.html)
  defining the type of each element in the array.

- `shape`:

  a vector of integers defining the shape of the array.

- `platform_config`:

  A [platform configuration](PlatformConfig.md) object

- `platform_config`:

  A [platform configuration](PlatformConfig.md) object

- `platform_config`:

  A [platform configuration](PlatformConfig.md) object

#### Returns

Returns the newly-created array stored at `key`.

------------------------------------------------------------------------

### Method [`print()`](https://rdrr.io/r/base/print.html)

Print-friendly representation of the object.

#### Usage

    SOMACollectionBase$print()

#### Returns

Invisibly returns `self`.

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    SOMACollectionBase$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
