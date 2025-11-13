# Base SOMA Context

Virtual R6 mapping class for SOMA context options. This class should be
used as the basis for platform-specific contexts as it checks
SOMA-specific context options

## See also

Derived classes: [`SOMATileDBContext`](SOMATileDBContext.md)

## Super classes

[`tiledbsoma::MappingBase`](MappingBase.md) -\>
[`tiledbsoma::ScalarMap`](ScalarMap.md) -\> `SOMAContextBase`

## Methods

### Public methods

- [`SOMAContextBase$new()`](#method-SOMAContextBase-new)

- [`SOMAContextBase$set()`](#method-SOMAContextBase-set)

- [`SOMAContextBase$clone()`](#method-SOMAContextBase-clone)

Inherited methods

- [`tiledbsoma::MappingBase$get()`](MappingBase.html#method-get)
- [`tiledbsoma::MappingBase$items()`](MappingBase.html#method-items)
- [`tiledbsoma::MappingBase$keys()`](MappingBase.html#method-keys)
- [`tiledbsoma::MappingBase$length()`](MappingBase.html#method-length)
- [`tiledbsoma::MappingBase$print()`](MappingBase.html#method-print)
- [`tiledbsoma::MappingBase$remove()`](MappingBase.html#method-remove)
- [`tiledbsoma::MappingBase$setv()`](MappingBase.html#method-setv)
- [`tiledbsoma::MappingBase$to_list()`](MappingBase.html#method-to_list)
- [`tiledbsoma::MappingBase$update()`](MappingBase.html#method-update)
- [`tiledbsoma::MappingBase$values()`](MappingBase.html#method-values)

------------------------------------------------------------------------

### Method `new()`

#### Usage

    SOMAContextBase$new(config = NULL)

#### Arguments

- `config`:

  ...

#### Returns

This is a **virtual** class and cannot be directly instantiated

------------------------------------------------------------------------

### Method `set()`

#### Usage

    SOMAContextBase$set(key, value)

#### Arguments

- `key`:

  Key to set

- `value`:

  Value to add for `key`, or `NULL` to remove the entry for `key`

#### Returns

\[chainable\] Invisibly returns `self` with `value` added for `key`

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    SOMAContextBase$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
