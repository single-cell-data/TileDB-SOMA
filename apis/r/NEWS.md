# 1.9.4

* This release contains an in-progress blockwise iterator, work on which will be completed in a subsequent release.

# 1.9.3

* This release contains a single Python-only performance improvment

# 1.9.2

* This release contains a single Python-only modification for an API backward-compatibility update

# 1.9.1

* This release contains a single Python-only modification for its build process

# 1.9.0

## Changes

* Add support for ingestion of `SeuratCommand` logs
* Add support for outgestion of `SeuratCommand` logs
* Add support for reading `*m` and `*p` layers from `SOMAExperimentAxisQuery`

# 1.7.0

## Changes

* Add support for registering arrays with their parent collection in `write_soma()`
* Add `write_soma.character()` to encode character vectors as SOMA arrays

# tiledbsoma 1.6.0

## Fixes

* `SOMADataFrame`, `SOMASparseNDArray`, and `SOMADenseNDArray`'s `write()` method now correctly leaves the array open in write mode

# tiledbsoma 1.5.0

## Changes

* Add support for I/O of R factors as enumerated types in `SOMADataFrame`
* Add support for writing `SummarizedExperiment` and `SingleCellExperiment` object to SOMAs
* Add support for bounding boxes for sparse arrays
* Add support for creating `SOMADataFrames` with `ordered()` columns
* Add support for I/O between SOMA and `SingleCellExperiment`
* Add support for updating `obs` and `var`
* Write string attrs as UTF-8 (Python compatibility)
* Optimize export of obsm/varm arrays to Seurat
* Add `axis_query()` method to `SOMAExperiment`
* Add `to_sparse_matrix()` method to `SOMAExperimentAxisQuery`
* Force exporting v3 assays with SeuratObject v5 installed

# tiledbsoma 1.4.0

This is the initial release of the tiledbsoma R package.

## Changes

* Added a `NEWS.md` file to track changes to the package
* `TileDBGroup` gains a `names` method to retrieve the names of group members
* Added `SOMAMeasurement` and `SOMAExperiment` classes
* spdl is now used for logging
* TileDB performance statistics can now be collected for analysis
* Added support for performing axis-based queries against a `SOMAExperiment` via the `SOMAExperimentAxisQuery` class
* `TileDBArray` class gained a `colnames()` method that returns the names of both dimensions and attributes
* Added internal helpers to centrally validate `coords` and `value_filter` arguments
* All R6 classes' `create()` method now return `self` rather than nothing
* Fixed calculating of relative paths when 1 of the URIs contains the `file://` prefix
* Added `PlatformConfig` and `SOMATileDBContext` classes to handle SOMA and TileDB configuration
* Add Seurat outgestors for `SOMAExperimentAxisQuery` objects
* Numeric coordinates passed to SOMADataFrame$read() are now automatically upcast to int64 when necessary
* Add ingestors to read data from `Seurat` objects
* Add methods for listing and accessing bundled datasets, which now includes a `SOMAExperiment` containing the pbmc_small dataset from the SeuratObject package
* New vignettes describing SOMA objects, reading data from them, and querying SOMA experiments
* Objects added to `SOMACollection`-based classes using the `add_new_*()` methods now pass through their parent context and platform config
* `SOMAExperimentAxisQuery` gained a `to_sparse_matrix()` method for retrieving data as a named sparse matrix
* `SOMAExperiment` gained `axis_query()` to construct a `SOMAExperimentAxisQuery` object
* Add SingleCellExperiment outgestor for `SOMAExperimentAxisQuery` objects
