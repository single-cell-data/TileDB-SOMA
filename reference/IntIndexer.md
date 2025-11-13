# The SOMA Re-Indexer

A re-indexer for unique integer indices

## Methods

### Public methods

- [`IntIndexer$new()`](#method-IntIndexer-new)

- [`IntIndexer$get_indexer()`](#method-IntIndexer-get_indexer)

- [`IntIndexer$clone()`](#method-IntIndexer-clone)

------------------------------------------------------------------------

### Method `new()`

Create a new re-indexer

#### Usage

    IntIndexer$new(data)

#### Arguments

- `data`:

  Integer keys used to build the index (hash) table

------------------------------------------------------------------------

### Method `get_indexer()`

Get the underlying indices for the target data

#### Usage

    IntIndexer$get_indexer(target, nomatch_na = FALSE)

#### Arguments

- `target`:

  Data to re-index

- `nomatch_na`:

  Set non-matches to `NA` instead of `-1`

#### Returns

A vector of 64-bit integers with `target` re-indexed

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    IntIndexer$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
(keys <- c(-10000, -100000, 200000, 5, 1, 7))
#> [1] -1e+04 -1e+05  2e+05  5e+00  1e+00  7e+00
(lookups <- unlist(replicate(n = 4L, c(-1L, 1:5), simplify = FALSE)))
#>  [1] -1  1  2  3  4  5 -1  1  2  3  4  5 -1  1  2  3  4  5 -1  1  2  3  4  5

indexer <- IntIndexer$new(keys)
indexer$get_indexer(lookups)
#> integer64
#>  [1] -1 4  -1 -1 -1 3  -1 4  -1 -1 -1 3  -1 4  -1 -1 -1 3  -1 4  -1 -1 -1 3 
indexer$get_indexer(lookups, nomatch_na = TRUE)
#> integer64
#>  [1] <NA> 4    <NA> <NA> <NA> 3    <NA> 4    <NA> <NA> <NA> 3    <NA> 4    <NA>
#> [16] <NA> <NA> 3    <NA> 4    <NA> <NA> <NA> 3   
```
