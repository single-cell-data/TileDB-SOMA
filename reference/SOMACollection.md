# SOMA Collection

Contains a key-value mapping where the keys are string names and the
values are any SOMA-defined foundational or composed type, including
`SOMACollection`, [`SOMADataFrame`](SOMADataFrame.md),
[`SOMADenseNDArray`](SOMADenseNDArray.md),
[`SOMASparseNDArray`](SOMASparseNDArray.md), or
[`SOMAExperiment`](SOMAExperiment.md) (lifecycle: maturing).

## Adding new objects to a collection

The `SOMACollection` class provides a number of type-specific methods
for adding new a object to the collection, such as
`add_new_sparse_ndarray()` and `add_new_dataframe()`. These methods will
create the new object and add it as member of the `SOMACollection`. The
new object will always inherit the parent context (see
[`SOMATileDBContext`](SOMATileDBContext.md)) and, by default, its
platform configuration (see [`PlatformConfig`](PlatformConfig.md)).
However, the user can override the default platform configuration by
passing a custom configuration to the `platform_config` argument.

## Carrara (TileDB v3) behavior

When working with Carrara URIs (`tiledb://workspace/teamspace/...`),
child objects created at a URI nested under a parent collection are
**automatically added** as members of the parent. This means:

- You do not need to call `add_new_collection()` after creating a child
  at a nested URI—the child is already a member.

- For backward compatibility, calling `add_new_collection()` on an
  already-registered child is a **no-op** and will not cause an error.

- The member name must match the relative URI segment (e.g., creating at
  `parent_uri/child` automatically adds the child with key `"child"`).

## Super classes

[`tiledbsoma::SOMAObject`](SOMAObject.md) -\>
[`tiledbsoma::SOMACollectionBase`](SOMACollectionBase.md) -\>
`SOMACollection`

## Methods

### Public methods

- [`SOMACollection$clone()`](#method-SOMACollection-clone)

Inherited methods

- [`tiledbsoma::SOMAObject$class()`](SOMAObject.html#method-class)
- [`tiledbsoma::SOMAObject$exists()`](SOMAObject.html#method-exists)
- [`tiledbsoma::SOMAObject$get_metadata()`](SOMAObject.html#method-get_metadata)
- [`tiledbsoma::SOMAObject$initialize()`](SOMAObject.html#method-initialize)
- [`tiledbsoma::SOMAObject$is_open()`](SOMAObject.html#method-is_open)
- [`tiledbsoma::SOMAObject$mode()`](SOMAObject.html#method-mode)
- [`tiledbsoma::SOMAObject$reopen()`](SOMAObject.html#method-reopen)
- [`tiledbsoma::SOMACollectionBase$add_new_collection()`](SOMACollectionBase.html#method-add_new_collection)
- [`tiledbsoma::SOMACollectionBase$add_new_dataframe()`](SOMACollectionBase.html#method-add_new_dataframe)
- [`tiledbsoma::SOMACollectionBase$add_new_dense_ndarray()`](SOMACollectionBase.html#method-add_new_dense_ndarray)
- [`tiledbsoma::SOMACollectionBase$add_new_sparse_ndarray()`](SOMACollectionBase.html#method-add_new_sparse_ndarray)
- [`tiledbsoma::SOMACollectionBase$close()`](SOMACollectionBase.html#method-close)
- [`tiledbsoma::SOMACollectionBase$create()`](SOMACollectionBase.html#method-create)
- [`tiledbsoma::SOMACollectionBase$get()`](SOMACollectionBase.html#method-get)
- [`tiledbsoma::SOMACollectionBase$length()`](SOMACollectionBase.html#method-length)
- [`tiledbsoma::SOMACollectionBase$names()`](SOMACollectionBase.html#method-names)
- [`tiledbsoma::SOMACollectionBase$open()`](SOMACollectionBase.html#method-open)
- [`tiledbsoma::SOMACollectionBase$print()`](SOMACollectionBase.html#method-print)
- [`tiledbsoma::SOMACollectionBase$remove()`](SOMACollectionBase.html#method-remove)
- [`tiledbsoma::SOMACollectionBase$set()`](SOMACollectionBase.html#method-set)
- [`tiledbsoma::SOMACollectionBase$set_metadata()`](SOMACollectionBase.html#method-set_metadata)

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    SOMACollection$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
uri <- withr::local_tempfile(pattern = "soma-collection")

(col <- SOMACollectionCreate(uri))
#> <SOMACollection>
#>   uri: /tmp/Rtmpfr8mYm/soma-collection2a7e66fa9b84
col$add_new_sparse_ndarray("sparse", arrow::float64(), shape = c(100L, 100L))
#> <SOMASparseNDArray>
#>   uri: /tmp/Rtmpfr8mYm/soma-collection2a7e66fa9b84/sparse
#>   dimensions: soma_dim_0, soma_dim_1 
#>   attributes: soma_data 
col$close()

(col <- SOMACollectionOpen(uri))
#> <SOMACollection>
#>   uri: /tmp/Rtmpfr8mYm/soma-collection2a7e66fa9b84
col$names()
#> [1] "sparse"
```
