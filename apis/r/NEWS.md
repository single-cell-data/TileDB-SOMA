# tiledbsoma 1.4.*

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
