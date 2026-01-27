# Ephemeral Collections

Ephemeral version of [`SOMACollection`](SOMACollection.md)`s`; ephemeral
collections are equivalent to [SOMA collections](SOMACollection.md) but
are stored in-memory instead of on-disk.

## Super classes

[`tiledbsoma::SOMAObject`](SOMAObject.md) -\>
[`tiledbsoma::SOMACollectionBase`](SOMACollectionBase.md) -\>
[`tiledbsoma::EphemeralCollectionBase`](EphemeralCollectionBase.md) -\>
`EphemeralCollection`

## Active bindings

- `soma_type`:

  The SOMA object type.

## Methods

### Public methods

- [`EphemeralCollection$clone()`](#method-EphemeralCollection-clone)

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

    EphemeralCollection$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
(col <- EphemeralCollection$new())
#> <EphemeralCollection>
#>   uri: ephemeral-collection:0x55cf36f51ef8
col$soma_type
#> [1] "SOMACollection"

dir <- withr::local_tempfile(pattern = "obs")
dir.create(dir, recursive = TRUE)

(obs <- load_dataset("soma-dataframe-pbmc3k-processed-obs", dir))
#> <SOMADataFrame>
#>   uri: /tmp/Rtmpfr8mYm/obs2a7e719aa1b1/soma-dataframe-pbmc3k-processed-obs
#>   dimensions: soma_joinid 
#>   attributes: orig.ident, nCount_RNA, nFeature_RNA, seurat_annotations, percent.mt, RNA_snn... 
col$set(obs, "obs")
col$names()
#> [1] "obs"
```
