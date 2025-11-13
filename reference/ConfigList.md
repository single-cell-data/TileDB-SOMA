# A Configuration List

An R6 mapping type for configuring various “parameters”. Essentially,
serves as a nested map where the inner map is a
[`ScalarMap`](ScalarMap.md):
`{<param>: `[`{<key>: <value>}`](ScalarMap.md)`}`

## Super class

[`tiledbsoma::MappingBase`](MappingBase.md) -\> `ConfigList`

## Methods

### Public methods

- [`ConfigList$get()`](#method-ConfigList-get)

- [`ConfigList$set()`](#method-ConfigList-set)

- [`ConfigList$setv()`](#method-ConfigList-setv)

- [`ConfigList$clone()`](#method-ConfigList-clone)

Inherited methods

- [`tiledbsoma::MappingBase$initialize()`](MappingBase.html#method-initialize)
- [`tiledbsoma::MappingBase$items()`](MappingBase.html#method-items)
- [`tiledbsoma::MappingBase$keys()`](MappingBase.html#method-keys)
- [`tiledbsoma::MappingBase$length()`](MappingBase.html#method-length)
- [`tiledbsoma::MappingBase$print()`](MappingBase.html#method-print)
- [`tiledbsoma::MappingBase$remove()`](MappingBase.html#method-remove)
- [`tiledbsoma::MappingBase$to_list()`](MappingBase.html#method-to_list)
- [`tiledbsoma::MappingBase$update()`](MappingBase.html#method-update)
- [`tiledbsoma::MappingBase$values()`](MappingBase.html#method-values)

------------------------------------------------------------------------

### Method [`get()`](https://rdrr.io/r/base/get.html)

#### Usage

    ConfigList$get(param, key = NULL, default = quote(expr = ))

#### Arguments

- `param`:

  Outer key or “parameter” to fetch

- `key`:

  Inner key to fetch; pass `NULL` to return the [map](ScalarMap.md) for
  `param`

- `default`:

  Default value to fetch if `key` is not found; defaults to `NULL`

#### Returns

The value of `key` for `param` in the map, or `default` if `key` is not
found

------------------------------------------------------------------------

### Method `set()`

#### Usage

    ConfigList$set(param, key, value)

#### Arguments

- `param`:

  Outer key or “parameter” to set

- `key`:

  Inner key to set

- `value`:

  Value to add for `key`, or `NULL` to remove the entry for `key`;
  optionally provide only `param` and `value` as a
  [`ScalarMap`](ScalarMap.md) to update `param` with the keys and values
  from `value`

#### Returns

\\chainable\\ Invisibly returns `self` with `value` added for `key` in
`param`

------------------------------------------------------------------------

### Method `setv()`

#### Usage

    ConfigList$setv(...)

#### Arguments

- `...`:

  Ignored

#### Returns

Nothing; `setv()` is disabled for `ConfigList` objects

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    ConfigList$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
(cfg <- ConfigList$new())
#> <ConfigList>
cfg$set("op1", "a", 1L)
cfg
#> <ConfigList>
#>   op1: <environment>
cfg$get("op1")
#> <ScalarMap>
#>   a: 1
```
