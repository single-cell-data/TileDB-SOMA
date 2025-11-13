# Platform Configuration

An R6 mapping type for configuring various “parameters” for multiple
“platforms”, essentially serves a multi-nested map where the inner map
is a [`ScalarMap`](ScalarMap.md) contained within a
[`ConfigList`](ConfigList.md) (middle map):
`{platform: {param: {key: value}}}`

## Super class

[`tiledbsoma::MappingBase`](MappingBase.md) -\> `PlatformConfig`

## Methods

### Public methods

- [`PlatformConfig$platforms()`](#method-PlatformConfig-platforms)

- [`PlatformConfig$params()`](#method-PlatformConfig-params)

- [`PlatformConfig$get()`](#method-PlatformConfig-get)

- [`PlatformConfig$get_params()`](#method-PlatformConfig-get_params)

- [`PlatformConfig$set()`](#method-PlatformConfig-set)

- [`PlatformConfig$setv()`](#method-PlatformConfig-setv)

- [`PlatformConfig$clone()`](#method-PlatformConfig-clone)

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

### Method `platforms()`

#### Usage

    PlatformConfig$platforms()

#### Returns

The names of the “platforms” (outer keys)

------------------------------------------------------------------------

### Method `params()`

#### Usage

    PlatformConfig$params(platform = NULL)

#### Arguments

- `platform`:

  The “platform” to pull parameter names (middle keys) for; pass `TRUE`
  to return all possible parameter names

#### Returns

The parameter names (middle keys) for `platform`

------------------------------------------------------------------------

### Method [`get()`](https://rdrr.io/r/base/get.html)

#### Usage

    PlatformConfig$get(
      platform,
      param = NULL,
      key = NULL,
      default = quote(expr = )
    )

#### Arguments

- `platform`:

  The name of the “platform” (outer key) to fetch

- `param`:

  The name of the “paramters” of `platform` to fetch; if `NULL`, returns
  the [configuration](ConfigList.md) for `platform`

- `key`:

  The “key” (inner key) for `param` in `platform` to fetch; if `NULL`
  and `param` is passed, returns the [map](ScalarMap.md) for `param` in
  `platform`

- `default`:

  Default value to fetch if `key` is not found; defaults to `null`

#### Returns

The value of `key` for `param` in `platform` in the map, or `default` if
`key` is not found

------------------------------------------------------------------------

### Method `get_params()`

#### Usage

    PlatformConfig$get_params(platform)

#### Arguments

- `platform`:

  The name of the “platform” (outer key) to fetch

#### Returns

The [`ConfigList`](ConfigList.md) for `platform`

------------------------------------------------------------------------

### Method `set()`

#### Usage

    PlatformConfig$set(platform, param, key, value)

#### Arguments

- `platform`:

  The name of the “platform” (outer key) to set

- `param`:

  Name of the “parameter” (middle key) in `platform` to set

- `key`:

  Inner key to set

- `value`:

  Value to add for `key`, or `NULL` to remove the entry for `key`;
  optionally provide only `platfomr`, `param`, and `value` as a
  [`ScalarMap`](ScalarMap.md) to update `param` for `platform` with the
  keys and values from `value`

#### Returns

\\chainable\\ Invisibly returns `self` with `value` added for `key` in
`param` for `platform`

------------------------------------------------------------------------

### Method `setv()`

#### Usage

    PlatformConfig$setv(...)

#### Arguments

- `...`:

  Ignored

#### Returns

Nothing; `setv()` is disabled for `PlatformConfig` objects

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    PlatformConfig$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
(cfg <- PlatformConfig$new())
#> <PlatformConfig>
cfg$set("plat1", "op1", "a", 1L)
cfg
#> <PlatformConfig>
#>   plat1: <environment>

cfg$get("plat1")
#> <ConfigList>
#>   op1: <environment>
cfg$get("plat1")$get("op1")
#> <ScalarMap>
#>   a: 1
```
