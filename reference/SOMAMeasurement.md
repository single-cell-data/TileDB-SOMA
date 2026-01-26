# SOMA Measurement

A `SOMAMeasurement` is a sub-element of a
[`SOMAExperiment`](SOMAExperiment.md), and is otherwise a specialized
[`SOMACollection`](SOMACollection.md) with pre-defined fields: `X`,
`var`, `obsm`/`varm`, and `obsp`/`varp` (see *Active Bindings* below for
details) (lifecycle: maturing).

## Adding new objects to a collection

The `SOMAMeasurement` class provides a number of type-specific methods
for adding new a object to the collection, such as
`add_new_sparse_ndarray()` and `add_new_dataframe()`. These methods will
create the new object and add it as member of the `SOMAMeasurement`. The
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
`SOMAMeasurement`

## Active bindings

- `var`:

  A [`SOMADataFrame`](SOMADataFrame.md) containing primary annotations
  on the variable axis, for variables in this measurement (i.e.,
  annotates columns of `X`). The contents of the `soma_joinid` column
  define the variable index domain, `var_id`. All variables for this
  measurement must be defined in this data frame.

- `X`:

  A [`SOMACollection`](SOMACollection.md) of
  [`SOMASparseNDArray`](SOMASparseNDArray.md)`s`, each contains measured
  feature values indexed by `[obsid, varid]`.

- `obsm`:

  A [`SOMACollection`](SOMACollection.md) of
  [`SOMADenseNDArray`](SOMADenseNDArray.md)`s` containing annotations on
  the observation axis. Each array is indexed by `obsid` and has the
  same shape as `obs`.

- `obsp`:

  A [`SOMACollection`](SOMACollection.md) of
  [`SOMASparseNDArray`](SOMASparseNDArray.md)`s` containing pairwise
  annotations on the observation axis and indexed with
  `[obsid_1, obsid_2]`.

- `varm`:

  A [`SOMACollection`](SOMACollection.md) of
  [`SOMADenseNDArray`](SOMADenseNDArray.md)`s` containing annotations on
  the variable axis. Each array is indexed by `varid` and has the same
  shape as `var`.

- `varp`:

  A [`SOMACollection`](SOMACollection.md) of
  [`SOMASparseNDArray`](SOMASparseNDArray.md)`s` containing pairwise
  annotations on the variable axis and indexed with
  `[varid_1, varid_2]`.

## Methods

### Public methods

- [`SOMAMeasurement$clone()`](#method-SOMAMeasurement-clone)

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

    SOMAMeasurement$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
uri <- withr::local_tempfile(pattern = "soma-measurement")
var <- data.frame(
  soma_joinid = bit64::seq.integer64(0L, 99L),
  var_id = paste0("feature_", seq_len(100L))
)
sch <- arrow::infer_schema(var)

(ms <- SOMAMeasurementCreate(uri))
#> <SOMAMeasurement>
#>   uri: /tmp/RtmpbAgXbM/soma-measurement284642c41e9c
sdf <- ms$add_new_dataframe(
  "var",
  sch,
  "soma_joinid",
  list(soma_joinid = c(0, 100))
)
sdf$write(arrow::as_arrow_table(var, schema = sch))
sdf$close()
ms$close()

(ms <- SOMAMeasurementOpen(uri))
#> <SOMAMeasurement>
#>   uri: /tmp/RtmpbAgXbM/soma-measurement284642c41e9c
ms$var
#> <SOMADataFrame>
#>   uri: file:///tmp/RtmpbAgXbM/soma-measurement284642c41e9c/var
#>   dimensions: soma_joinid 
#>   attributes: var_id 
```
