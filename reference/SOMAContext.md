# SOMA Context

Context map for TileDB-SOMA objects

## Active bindings

- `handle`:

  External pointer to the C++ interface

## Methods

### Public methods

- [`SOMAContext$new()`](#method-SOMAContext-new)

- [`SOMAContext$get_config()`](#method-SOMAContext-get_config)

- [`SOMAContext$get_data_protocol()`](#method-SOMAContext-get_data_protocol)

- [`SOMAContext$is_tiledbv2()`](#method-SOMAContext-is_tiledbv2)

- [`SOMAContext$is_tiledbv3()`](#method-SOMAContext-is_tiledbv3)

- [`SOMAContext$clone()`](#method-SOMAContext-clone)

------------------------------------------------------------------------

### Method `new()`

#### Usage

    SOMAContext$new(config = NULL)

#### Arguments

- `config`:

  ...

#### Returns

An instantiated `SOMATileDBContext` object

------------------------------------------------------------------------

### Method `get_config()`

#### Usage

    SOMAContext$get_config()

#### Returns

A character vector with the current config options set on the context.

------------------------------------------------------------------------

### Method `get_data_protocol()`

#### Usage

    SOMAContext$get_data_protocol(uri)

#### Arguments

- `uri`:

  A URI for a SOMA object

#### Returns

The data protocol to use for the URI.

------------------------------------------------------------------------

### Method `is_tiledbv2()`

#### Usage

    SOMAContext$is_tiledbv2(uri)

#### Arguments

- `uri`:

  A URI for a SOMA object

#### Returns

TRUE if the URI will use `tiledbv2` semantics.

------------------------------------------------------------------------

### Method `is_tiledbv3()`

#### Usage

    SOMAContext$is_tiledbv3(uri)

#### Arguments

- `uri`:

  A URI for a SOMA object

#### Returns

TRUE if the URI will use `tiledbv3` semantics.

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    SOMAContext$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
