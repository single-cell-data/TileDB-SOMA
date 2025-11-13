# A Mapping Type for Scalars

An R6 mapping type that is limited to scalar atomic vector types only;
can optionally be limited further to a specific atomic vector type (eg.
“`logical`”).

## Super class

[`tiledbsoma::MappingBase`](MappingBase.md) -\> `ScalarMap`

## Active bindings

- `type`:

  The type that this `ScalarMap` is limited to

## Methods

### Public methods

- [`ScalarMap$new()`](#method-ScalarMap-new)

- [`ScalarMap$set()`](#method-ScalarMap-set)

- [`ScalarMap$clone()`](#method-ScalarMap-clone)

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

    ScalarMap$new(type = "any")

#### Arguments

- `type`:

  Limit the `ScalarMap` to a preset type; choose from:

  - “`any`”

  - “`numeric`”

  - “`integer`”

  - “`character`”

  - “`logical`”

#### Returns

An instantiated `ScalarMap` object with the type set to `type`

------------------------------------------------------------------------

### Method `set()`

#### Usage

    ScalarMap$set(key, value)

#### Arguments

- `key`:

  Key to set

- `value`:

  Value to add for `key`, or `NULL` to remove the entry for `key`

#### Returns

\[chainable\] Invisibly returns `self` with `value` added as `key`

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    ScalarMap$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
(map <- ScalarMap$new())
#> <ScalarMap>
map$set("a", 1L)
map
#> <ScalarMap>
#>   a: 1

map$get("a")
#> [1] 1
map$get("b", default = NULL)
#> NULL
```
