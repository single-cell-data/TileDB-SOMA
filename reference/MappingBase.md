# R6 Base Mapping Type

Virtual base mapping type for R6 objects; defines internal data
structure (`private$.data`) as a named list along with behavior methods
for getting (`self$get()`) and setting (`self$set()`) items in the map

## Usage

``` r
# S3 method for class 'MappingBase'
x[[i, ..., default = NULL]]

# S3 method for class 'MappingBase'
x[[i, ...]] <- value

# S3 method for class 'MappingBase'
as.list(x, ...)

# S3 method for class 'MappingBase'
length(x)

# S3 method for class 'MappingBase'
names(x)
```

## Arguments

- x:

  A mapping object

- i:

  A key to fetch or set; see `$get()` or `$set()` methods below

- ...:

  Ignored

- default:

  Default value to fetch if `i` is not found; defaults to `NULL`

- value:

  Value to add for `i`, or `NULL` to remove the entry for `i`

## Value

`[[`: The value of `i` in the map, or `default` if `i` is not found

`[[<-`: `x` with `value` added as `i`

`as.list`: The map as a list

`length`: The number of items in the map

`names`: The keys of the map

## See also

Derived classes: [`ConfigList`](ConfigList.md),
[`PlatformConfig`](PlatformConfig.md), [`ScalarMap`](ScalarMap.md),
[`TileDBCreateOptions`](TileDBCreateOptions.md)

## Methods

### Public methods

- [`MappingBase$new()`](#method-MappingBase-new)

- [`MappingBase$keys()`](#method-MappingBase-keys)

- [`MappingBase$values()`](#method-MappingBase-values)

- [`MappingBase$items()`](#method-MappingBase-items)

- [`MappingBase$get()`](#method-MappingBase-get)

- [`MappingBase$set()`](#method-MappingBase-set)

- [`MappingBase$setv()`](#method-MappingBase-setv)

- [`MappingBase$remove()`](#method-MappingBase-remove)

- [`MappingBase$update()`](#method-MappingBase-update)

- [`MappingBase$length()`](#method-MappingBase-length)

- [`MappingBase$to_list()`](#method-MappingBase-to_list)

- [`MappingBase$print()`](#method-MappingBase-print)

- [`MappingBase$clone()`](#method-MappingBase-clone)

------------------------------------------------------------------------

### Method `new()`

#### Usage

    MappingBase$new(...)

#### Arguments

- `...`:

  Ignored

#### Returns

This is a **virtual** class and cannot be directly instantiated

------------------------------------------------------------------------

### Method `keys()`

#### Usage

    MappingBase$keys()

#### Returns

The keys of the map

------------------------------------------------------------------------

### Method `values()`

#### Usage

    MappingBase$values()

#### Returns

A `list` containing the map values

------------------------------------------------------------------------

### Method `items()`

#### Usage

    MappingBase$items()

#### Returns

Return the items of the map as a list

------------------------------------------------------------------------

### Method [`get()`](https://rdrr.io/r/base/get.html)

#### Usage

    MappingBase$get(key, default = quote(expr = ))

#### Arguments

- `key`:

  Key to fetch

- `default`:

  Default value to fetch if `key` is not found; defaults to `NULL`

#### Returns

The value of `key` in the map, or `default` if `key` is not found

------------------------------------------------------------------------

### Method `set()`

#### Usage

    MappingBase$set(key, value)

#### Arguments

- `key`:

  Key to set

- `value`:

  Value to add for `key`, or `NULL` to remove the entry for `key`

#### Returns

\[chainable\] Invisibly returns `self` with `value` added as `key`

------------------------------------------------------------------------

### Method `setv()`

#### Usage

    MappingBase$setv(...)

#### Arguments

- `...`:

  Named arguments to add to `self`

#### Returns

\[chainable\] Invisibly returns `self` with the values of `...` added to
the map

------------------------------------------------------------------------

### Method [`remove()`](https://rdrr.io/r/base/rm.html)

#### Usage

    MappingBase$remove(key)

#### Arguments

- `key`:

  Key to remove

#### Returns

\[chainable\] Invisibly returns `self` with `key` removed from the map

------------------------------------------------------------------------

### Method [`update()`](https://rdrr.io/r/stats/update.html)

#### Usage

    MappingBase$update(map)

#### Arguments

- `map`:

  A mapping type to update the current map with

#### Returns

\[chainable\] Invisibly returns `self` with the value of `map`

------------------------------------------------------------------------

### Method [`length()`](https://rdrr.io/r/base/length.html)

#### Usage

    MappingBase$length()

#### Returns

The number of items in the map

------------------------------------------------------------------------

### Method `to_list()`

#### Usage

    MappingBase$to_list()

#### Returns

The map as a list

------------------------------------------------------------------------

### Method [`print()`](https://rdrr.io/r/base/print.html)

#### Usage

    MappingBase$print()

#### Returns

\[chainable\] Prints information about the map to `stdout` and invisibly
returns `self`

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    MappingBase$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
