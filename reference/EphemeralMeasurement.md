# Ephemeral SOMA Measurement

Ephemeral version of [`SOMAMeasurement`](SOMAMeasurement.md)`s`;
ephemeral measurements are equivalent to [SOMA
measurements](SOMAMeasurement.md) but are stored in-memory instead of
on-disk.

## Super classes

[`tiledbsoma::SOMAObject`](SOMAObject.md) -\>
[`tiledbsoma::SOMACollectionBase`](SOMACollectionBase.md) -\>
[`tiledbsoma::EphemeralCollectionBase`](EphemeralCollectionBase.md) -\>
`EphemeralMeasurement`

## Active bindings

- `var`:

  A [`SOMADataFrame`](SOMADataFrame.md) containing primary annotations
  on the variable axis, for variables in this measurement (i.e.,
  annotates columns of `X`). The contents of the `soma_joinid` column
  define the variable index domain, `var_id`. All variables for this
  measurement must be defined in this data frame.

- `X`:

  A [`SOMACollection`](SOMACollection.md) of
  [`SOMASparseNDArray`](SOMASparseNDArray.md)s; each contain measured
  feature values indexed by `[obsid, varid]`.

- `obsm`:

  A [`SOMACollection`](SOMACollection.md) of
  [`SOMADenseNDArray`](SOMADenseNDArray.md)s containing annotations on
  the observation axis. Each array is indexed by `obsid` and has the
  same shape as `obs`.

- `obsp`:

  A [`SOMACollection`](SOMACollection.md) of
  [`SOMASparseNDArray`](SOMASparseNDArray.md)s containing pairwise
  annotations on the observation axis and indexed with
  `[obsid_1, obsid_2]`.

- `varm`:

  A [`SOMACollection`](SOMACollection.md) of
  [`SOMADenseNDArray`](SOMADenseNDArray.md)s containing annotations on
  the variable axis. Each array is indexed by `varid` and has the same
  shape as `var`.

- `varp`:

  A [`SOMACollection`](SOMACollection.md) of
  [`SOMASparseNDArray`](SOMASparseNDArray.md)s containing pairwise
  annotations on the variable axis and indexed with
  `[varid_1, varid_2]`.

- `soma_type`:

  The SOMA object type.

## Methods

### Public methods

- [`EphemeralMeasurement$clone()`](#method-EphemeralMeasurement-clone)

Inherited methods

- [`tiledbsoma::SOMAObject$class()`](SOMAObject.html#method-class)
- [`tiledbsoma::SOMAObject$is_open()`](SOMAObject.html#method-is_open)
- [`tiledbsoma::SOMAObject$mode()`](SOMAObject.html#method-mode)
- [`tiledbsoma::SOMAObject$reopen()`](SOMAObject.html#method-reopen)
- [`tiledbsoma::EphemeralCollectionBase$add_new_collection()`](EphemeralCollectionBase.html#method-add_new_collection)
- [`tiledbsoma::EphemeralCollectionBase$add_new_dataframe()`](EphemeralCollectionBase.html#method-add_new_dataframe)
- [`tiledbsoma::EphemeralCollectionBase$add_new_dense_ndarray()`](EphemeralCollectionBase.html#method-add_new_dense_ndarray)
- [`tiledbsoma::EphemeralCollectionBase$add_new_sparse_ndarray()`](EphemeralCollectionBase.html#method-add_new_sparse_ndarray)
- [`tiledbsoma::EphemeralCollectionBase$close()`](EphemeralCollectionBase.html#method-close)
- [`tiledbsoma::EphemeralCollectionBase$create()`](EphemeralCollectionBase.html#method-create)
- [`tiledbsoma::EphemeralCollectionBase$exists()`](EphemeralCollectionBase.html#method-exists)
- [`tiledbsoma::EphemeralCollectionBase$get()`](EphemeralCollectionBase.html#method-get)
- [`tiledbsoma::EphemeralCollectionBase$get_metadata()`](EphemeralCollectionBase.html#method-get_metadata)
- [`tiledbsoma::EphemeralCollectionBase$get_tiledb_config()`](EphemeralCollectionBase.html#method-get_tiledb_config)
- [`tiledbsoma::EphemeralCollectionBase$initialize()`](EphemeralCollectionBase.html#method-initialize)
- [`tiledbsoma::EphemeralCollectionBase$length()`](EphemeralCollectionBase.html#method-length)
- [`tiledbsoma::EphemeralCollectionBase$names()`](EphemeralCollectionBase.html#method-names)
- [`tiledbsoma::EphemeralCollectionBase$open()`](EphemeralCollectionBase.html#method-open)
- [`tiledbsoma::EphemeralCollectionBase$print()`](EphemeralCollectionBase.html#method-print)
- [`tiledbsoma::EphemeralCollectionBase$remove()`](EphemeralCollectionBase.html#method-remove)
- [`tiledbsoma::EphemeralCollectionBase$set()`](EphemeralCollectionBase.html#method-set)
- [`tiledbsoma::EphemeralCollectionBase$set_metadata()`](EphemeralCollectionBase.html#method-set_metadata)

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    EphemeralMeasurement$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
(ms <- EphemeralMeasurement$new())
#> <EphemeralMeasurement>
#>   uri: ephemeral-collection:0x558432db9cc8
ms$soma_type
#> [1] "SOMAMeasurement"

ms$set(EphemeralCollection$new(), "X")
ms$X
#> <EphemeralCollection>
#>   uri: ephemeral-collection:0x558432e73cc0
```
