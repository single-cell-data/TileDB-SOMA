# TileDB-SOMA Python Changelog

All notable changes to the Python TileDB-SOMA project will be documented in this file (related: [TileDB-SOMA R API changelog](../r/NEWS.md)).

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [Unreleased]

### Added

### Changed

- \[[#4126](https://github.com/single-cell-data/TileDB-SOMA/pull/4126)\] [python] At package import time, validate that the expected TileDB version is installed and used. Raises a RuntimeError exception if the condition is not met. This is an attempt to better warn users who have corrupted conda installations.
- \[[#4137](https://github.com/single-cell-data/TileDB-SOMA/pull/4137)\] [python][BREAKING] Add user-specified obs/var index column to `ExperimentAxisQuery.to_anndata`. This change will also set the index columns based upon metadata hints, following the conventions of `tiledbsoma.io`. See the docstrings for more details. Prior to this change, the index dtype would be set to `string` in all cases - this change removes the forced cast, and leaves the column type as-is.
- Update TileDB version to https://github.com/single-cell-data/TileDB-SOMA/pull/4177

### Deprecated

### Removed

### Fixed

### Security

## [Release 1.17.0]

Release 1.17.0 upgrades from TileDB 2.28.0 to TileDB 2.28.1.


### Changed

- \[[#4177](https://github.com/single-cell-data/TileDB-SOMA/pull/4177)\] Update [TileDB core to 2.28.1](https://github.com/TileDB-Inc/TileDB/releases/tag/2.28.1).


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

## [Release prior to 1.17.0]

TileDB-SOMA Python releases prior to 1.17.0 are documented in the [TileDB-SOMA Github Releases](https://github.com/single-cell-data/TileDB-SOMA/releases).
