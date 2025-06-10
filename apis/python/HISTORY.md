# TileDB-SOMA Python Changelog

All notable changes to the Python TileDB-SOMA project will be documented in this file (related: [TileDB-SOMA R API changelog](../r/NEWS.md)).

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [Unreleased]

### Added

### Changed

- \[[#3983](https://github.com/single-cell-data/TileDB-SOMA/pull/3983)\] [python] Multiple writes of pre-sorted data may now be written to a single fragment using TileDB global order writes. Enable this performance optimization by setting the platform_config parameter `sort_coords` to `False` in the call to write. Will raise an error if data is not written in global sort order.
- \[[#4086](https://github.com/single-cell-data/TileDB-SOMA/pull/4086)\] [python] Add new parameter `allow_duplicate_obs_ids` to the `tiledbsoma.io` functions `register_anndatas` and `register_h5ads`.  When `False` (default), a error will be raised if there are any duplicate `obs` IDs in the provided SOMA Experiment or AnnData objects. Set the parameter to `True` for legacy behavior. ID handling on the `var` axis is unchanged.
- \[[#4108](https://github.com/single-cell-data/TileDB-SOMA/pull/4108)\] [python] improve performance of `tiledbsoma.io.from_anndata` and `from_h5ad` when appending groups of AnnData known to have no duplicate obs axis IDs.

- \[[#4106](https://github.com/single-cell-data/TileDB-SOMA/pull/4106)\] [python][BREAKING] The `SOMAObject.reopen` method now modifies the orginal `SOMAObject` in place (flushes data to disk and reopens with the requested timestamp and mode) and returns a reference to itself instead of flushing data to disk and opening a new object.

### Deprecated

- \[[#4081](https://github.com/single-cell-data/TileDB-SOMA/pull/4081)\] [python] the `tiledbsoma.io` functions `append_obs`, `append_var` and `append_X` are deprecated and will be removed in a future release. It is recommended to use tiledbsoma.io.from_anndata (with a registration map from tiledbsoma.io.register_anndatas or tiledbsoma.io.register_h5ads) for appending new, complete AnnData objects to an Experiment.
- \[[#4082](https://github.com/single-cell-data/TileDB-SOMA/pull/4082)\] [python] `tiledbsoma.io.create_from_matrix` is deprecated and will be removed in a future release. To add a new matrix as a layer within an existing SOMA Experiment (e.g., to X, obsm, varm), please use the more specific functions tiledbsoma.io.add_X_layer or tiledbsoma.io.add_matrix_to_collection. If you need to create a standalone SOMA NDArray outside of a pre-defined Experiment structure, please use the direct SOMA API constructors, such as tiledbsoma.SparseNDArray.create.
- \[[#4083](https://github.com/single-cell-data/TileDB-SOMA/pull/4083)\] [python] "resume" mode in tiledbsoma.io ingestion methods is deprecated and will be removed i a future release. This includes from_anndata, from_h5ad and related ingest functions. The recommended and safest approach for recovering from a failed ingestion is to delete the partially written SOMA Experiment and restart the ingestion process from the original input files or a known-good backup.

### Removed

### Fixed

- \[[#4071](https://github.com/single-cell-data/TileDB-SOMA/pull/4071)\] [python] A `tiledb_timestamp` with value of zero is now equivalent to an unspecified timestamp (or `None`), and will be a synonym for "current time". Prior to this fix, a zero-valued timestamp would generate errors or unpredictable results.
- \[[#4103](https://github.com/single-cell-data/TileDB-SOMA/pull/4103)\] [python] Do not attempt to resize empty measurements in tiledbsoma.io.prepare_experiment. This fixes a bug where calling `prepare_experiment` would fail if the Experiment contained any empty Measurements.

### Security

## [Release 1.17.0]

The primary change in 1.17.0 is the upgrade to TileDB 2.28.

### Added

- \[[#3740](https://github.com/single-cell-data/TileDB-SOMA/pull/3740)\] [python] Add experimental Dask-backed `to_anndata` functionality to `SparseNDArray` and `ExperimentAxisQuery`.

### Changed

- \[\[[#4057](https://github.com/single-cell-data/TileDB-SOMA/pull/4057)\] [c++] Update [TileDB core to 2.28.0](https://github.com/TileDB-Inc/TileDB/blob/main/HISTORY.md#tiledb-v2280-release-notes)
- \[[#4023](https://github.com/single-cell-data/TileDB-SOMA/pull/4023)\] [c++] Use nanoarrow ArrowSchemaSetTypeDateTime for datetime values. Dictionary type with timestamp value type will raise error on read.

### Fixed

- \[[#4040](https://github.com/single-cell-data/TileDB-SOMA/pull/4040)\] [python] suppress insignificant overflow warning from numpy.
- \[[#4050](https://github.com/single-cell-data/TileDB-SOMA/pull/4050)\] DataFrame `count` and SparseNDArray `nnz` fix - report correct number of cells in array in the case where a delete query had been previously applied.
- \[[#4055](https://github.com/single-cell-data/TileDB-SOMA/pull/4055)\] Various `open()` code paths failed to check the SOMA encoding version number, and would fail with cryptic errors.
- \[[#4066](https://github.com/single-cell-data/TileDB-SOMA/pull/4066)\] Fix various memory leaks related to releasing Arrow structures when transfering ownership between C++ and Python and vise versa.
- \[[#4031](https://github.com/single-cell-data/TileDB-SOMA/pull/4031)\] [python] Storage paths generated from collection keys are now URL-escaped if they contain characters outside the safe set (`a-zA-Z0-9-_.()^!@+={}~'`). Additionally, the special names `..` and `.` are now prohibited.

## [Release prior to 1.17.0]

TileDB-SOMA Python releases prior to 1.17.0 are documented in the [TileDB-SOMA Github Releases](https://github.com/single-cell-data/TileDB-SOMA/releases).
