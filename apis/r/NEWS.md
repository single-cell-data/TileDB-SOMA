# Unreleased

## Added

- Use a global memory budget for read operations instead of a per column memory budget. The global memory budget allocates splits the budget per column depending on the type and characteristics of each column. Global memory budget is disabled by default under a feature flag and can be enabled by setting `soma.read.use_memory_pool`. ([#4299](https://github.com/single-cell-data/TileDB-SOMA/pull/4299))

## Changed

- `ManagedQuery` reuses the same buffers for each incomplete read and allocates dedicated buffers when converting to Arrow. ([#4299](https://github.com/single-cell-data/TileDB-SOMA/pull/4299))

## Deprecated

- The function `soma_context` is deprecated. Use class `SOMAContext` instead. ([#4355](https://github.com/single-cell-data/TileDB-SOMA/pull/4355))
- The parameter `tiledbsoma_ctx` is deprecated in all functions/methods that use it. Use the parameter `context` instead. ([#4355](https://github.com/single-cell-data/TileDB-SOMA/pull/4355))
- `SOMATileDBContext` is deprecated. Use class `SOMAContext` instead. ([#4355](https://github.com/single-cell-data/TileDB-SOMA/pull/4355))

## Removed

## Fixed

- The SOMA Context is only cached as an environment variable when the function `soma_context` is called directly. ([#4355](https://github.com/single-cell-data/TileDB-SOMA/pull/4355))
- `SOMATileDBContext` no longer replaces `sm.mem.reader.sparse_global_order.ratio_array_data` when set in the input config. ([#4355](https://github.com/single-cell-data/TileDB-SOMA/pull/4355))

## Security

# tiledbsoma 2.2.0

## Changed

- Update TileDB core version to [2.30.0](https://github.com/TileDB-Inc/TileDB/releases/tag/2.30.0). ([#4324](https://github.com/single-cell-data/TileDB-SOMA/pull/4324))

# tiledbsoma 2.1.2

**Action required:** Users who ingested BPCells data with version 2.1.0 or 2.1.1 must re-ingest with this version to ensure data correctness.

## Fixed

- Addressed bug in BPCells ingestion, where data ingested in an iterated manner would only be written at the top of the array. Ingestion now proceeds down the array for each iteration. ([#4317](https://github.com/single-cell-data/TileDB-SOMA/pull/4317))

# tiledbsoma 2.1.1

## Changed

- Update TileDB core version to [2.29.2](https://github.com/TileDB-Inc/TileDB/releases/tag/2.29.2). ([#4309](https://github.com/single-cell-data/TileDB-SOMA/pull/4309))

# tiledbsoma 2.1.0

This release adds support for ingestion of BPCells-backed `Seurat` objects in `write_soma()`, adds warnings for new deprecations, and updates the TileDB core version to 2.29.1.

## Added

- Add support for ingestion of BPCells-backed `Seurat` objects in `write_soma()` ([#4273](https://github.com/single-cell-data/TileDB-SOMA/pull/4273))

## Changed

- Update TileDB core version to [2.29.1](https://github.com/TileDB-Inc/TileDB/releases/tag/2.29.1). ([#4284](https://github.com/single-cell-data/TileDB-SOMA/pull/4284))

## Deprecated

- Mark `SOMANDArrayBase$set_data_type()` for deprecation. ([#4279](https://github.com/single-cell-data/TileDB-SOMA/pull/4279))
- Mark setting `domain` to `NULL` in `SOMADataFrame$create()` for deprecation. ([#4274](https://github.com/single-cell-data/TileDB-SOMA/pull/4274))

# tiledbsoma 2.0.0

This release is the first TileDB-SOMA release that follows our new versioning policy (see the [developer docs](../../dev_docs/POLICIES.md)). It introduces the new `DELETE` mode and deprecates removing elements from a `Collection` in `WRITE` mode. It also updates the core TileDB version to 2.28.1.

## Added

- Add delete mode specified by `mode="DELETE"`. ([#4125](https://github.com/single-cell-data/tiledb-soma/pull/4125))

## Changed

- Update [TileDB core to 2.28.1](https://github.com/TileDB-Inc/TileDB/releases/tag/2.28.1). ([#4077](https://github.com/single-cell-data/TileDB-SOMA/pull/4177))
- Remove RcppSpdlog and spld as dependencies. Logs are no longer forwarded through R output.

## Deprecated

- Deprecate removing elements from a collection in write mode. In the future, all new removals will need to be done in delete mode. ([#4125](https://github.com/single-cell-data/tiledb-soma/pull/4125))

# tiledbsoma 1.18.0

## Added

- New `SOMAObject` base class to serve as root for `SOMAArrayBase` and `SOMACollectionBase` ([#3771](https://github.com/single-cell-data/TileDB-SOMA/pull/3771))
- New field `SOMACollectionBase$members` to get a list with the members of a collection ([#3771](https://github.com/single-cell-data/TileDB-SOMA/pull/3771))

## Changed

- Handshake `internal_use_only = "allowed_use"` for `$new()`, `$open()`, and `$create()` has been replaced with environment scoping; use of factory functions for opening and creation is now mandatory ([#3771](https://github.com/single-cell-data/TileDB-SOMA/pull/3771))
- `TileDBObject`, `TileDBArray`, and `TileDBGroup` classes have been removed ([#3771](https://github.com/single-cell-data/TileDB-SOMA/pull/3771))
- `SOMACollection$to_list()` and `SOMACollection$to_data_frame()` have been removed as they were unused public internal methods inherited from `TileDBGroup` ([#3771](https://github.com/single-cell-data/TileDB-SOMA/pull/3771))

## Removed

- `SOMAArray$object` has been removed ([#3771](https://github.com/single-cell-data/TileDB-SOMA/pull/3771))

# tiledbsoma 1.17.0

## Removed

- `TileDBURI` class has been removed ([#3644](https://github.com/single-cell-data/TileDB-SOMA/pull/3644))
- arrow \<> TileDB-R helpers have been removed ([#3644](https://github.com/single-cell-data/TileDB-SOMA/pull/3644))
  - `arrow_field_from_tiledb_dim()`
  - `arrow_field_from_tiledb_attr()`
  - `arrow_schema_from_tiledb_schema()`

## Added

- New function `get_tiledb_object_type()` to replace `tiledb::tiledb_object_type()` ([#3644](https://github.com/single-cell-data/TileDB-SOMA/pull/3644))
- New function `get_tiledb_version()` to replace `tiledb::tiledb_version()` ([#3644](https://github.com/single-cell-data/TileDB-SOMA/pull/3644))
- New method `SOMAArrayBase$is_sparse()` to replace `tiledb::is.sparse()` ([#3644](https://github.com/single-cell-data/TileDB-SOMA/pull/3644))
- New method `SOMAArrayBase$allows_duplicates()` to replace `tiledb::allows_dups()` ([#3644](https://github.com/single-cell-data/TileDB-SOMA/pull/3644))
- New method `SOMADataFrame$levels()` to replace tiledb-r enum accessors ([#3644](https://github.com/single-cell-data/TileDB-SOMA/pull/3644))

## Changed

- Update [TileDB core to 2.28.0](https://github.com/TileDB-Inc/TileDB/blob/main/HISTORY.md#tiledb-v2280-release-notes) ([#4057](https://github.com/single-cell-data/TileDB-SOMA/pull/4057))
- `TileDBArray$attributes()` has been promoted to `SOMAArrayBase$attributes()` and returns a named list instead of an external pointer ([#3644](https://github.com/single-cell-data/TileDB-SOMA/pull/3644))

## Fixed

- Fix C++20 flag to be a configuration option instead of hard-coded ([#4051](https://github.com/single-cell-data/TileDB-SOMA/pull/4051))

# tiledbsoma 1.16.0

- Encode string metadata as `TILEDB_STRING_UTF8` instead of `TILEDB_STRING_ASCII`
- Use S3 method dispatch on `integer64` instead of directly calling the S3 methods
- \[c++\] Replace `SOMAArray` read and write calls with `ManagedQuery` [#3678](https://github.com/single-cell-data/TileDB-SOMA/pull/3678)
- Remove `used_shape`, which was deprecated in 1.15 with scheduled removal in 1.16 [#3723](https://github.com/single-cell-data/TileDB-SOMA/pull/3723)

# tiledbsoma 1.15.7

- This release contains Python-only bugfixes

# tiledbsoma 1.15.6

- This release contains Python-only bugfixes

# tiledbsoma 1.15.5

- This release contains a Python-only bugfix

# tiledbsoma 1.15.4

- This release contains Python-only updates for spatial transcriptomics

# tiledbsoma 1.15.3

- This release contains Python-only bugfixes

# tiledbsoma 1.15.2

- This release contains a Python-only bugfix

# tiledbsoma 1.15.1

- Encode string metadata as `TILEDB_STRING_UTF8` instead of `TILEDB_STRING_ASCII` [#3469](https://github.com/single-cell-data/TileDB-SOMA/pull/3469)

# tiledbsoma 1.15.0

## Changes

- Remove unused `fragment_count` accessor [#3054](https://github.com/single-cell-data/TileDB-SOMA/pull/3054)
- Implement missing `domain` argument to `SOMADataFrame` `create` [#3032](https://github.com/single-cell-data/TileDB-SOMA/pull/3032)
- Use `libtiledbsoma` for R schema evolution [#3100](https://github.com/single-cell-data/TileDB-SOMA/pull/3100)
- Push `attrnames` down to C++ [#3121](https://github.com/single-cell-data/TileDB-SOMA/pull/3121)
- Push `schema` accessor down to `libtiledbsoma` [#3079](https://github.com/single-cell-data/TileDB-SOMA/pull/3079)
- Handle `numeric` coords properly when reading arrays
- Remove two more `tiledb::schema` callsites [#3160](https://github.com/single-cell-data/TileDB-SOMA/pull/3160)
- Add new Arrow-to-R type mapper
- Add transitional/non-exported `parse_query_condition_new` [#3162](https://github.com/single-cell-data/TileDB-SOMA/pull/3162)
- Apply new `parse_query_condition` [#3174](https://github.com/single-cell-data/TileDB-SOMA/pull/3174)
- Apply new `non_empty_domain` [#3176](https://github.com/single-cell-data/TileDB-SOMA/pull/3176)
- Support for dense current domain with core 2.27 [#3180](https://github.com/single-cell-data/TileDB-SOMA/pull/3180)
- Fix `is_named_list` bug for half-named lists [#3183](https://github.com/single-cell-data/TileDB-SOMA/pull/3183)
- Expose block/random writer for sparse arrays [#3204](https://github.com/single-cell-data/TileDB-SOMA/pull/3204)
- Min-sizing for dataframes/arrays with new shape feature [#3208](https://github.com/single-cell-data/TileDB-SOMA/pull/3208)
- Proper prefixing for shape-related methods [#3237](https://github.com/single-cell-data/TileDB-SOMA/pull/3237)
- Bindings for `upgrade_domain` [#3238](https://github.com/single-cell-data/TileDB-SOMA/pull/3238)
- Apply `styler::style_pkg()` [#3239](https://github.com/single-cell-data/TileDB-SOMA/pull/3239)
- Plumb old-style `SOMATileDBContext` into new-style `soma_context()` [#3252](https://github.com/single-cell-data/TileDB-SOMA/pull/3252)
- Fixes for dense arrays and yet-to-be-released core 2.27 [#3270](https://github.com/single-cell-data/TileDB-SOMA/pull/3270)
- More fixes for unit-test cases with dense + core 2.27 [#3280](https://github.com/single-cell-data/TileDB-SOMA/pull/3280)
- Add initial support for ragged array writing for Seurat v5 by @mojaveazure in [#2523](https://github.com/single-cell-data/TileDB-SOMA/pull/2523)
- Remove 2.27-related feature flag by @johnkerl in [#3368](https://github.com/single-cell-data/TileDB-SOMA/pull/3368)
- Revert #3300 by @johnkerl in [#3358](https://github.com/single-cell-data/TileDB-SOMA/pull/3358)
- Enforce dataframe domain lower bound == 0 by @johnkerl in [#3300](https://github.com/single-cell-data/TileDB-SOMA/pull/3300)
- Data refresh for new shapes by @johnkerl in [#3303](https://github.com/single-cell-data/TileDB-SOMA/pull/3303)
- Move beyond the new-shape feature flag by @johnkerl in [#3301](https://github.com/single-cell-data/TileDB-SOMA/pull/3301)
- Note on dev installs and `spdlog` failures by @johnkerl in [#3324](https://github.com/single-cell-data/TileDB-SOMA/pull/3324)
- Remove tiledb-r specific install by @mojaveazure in [#3319](https://github.com/single-cell-data/TileDB-SOMA/pull/3319)
- Try to unbreak `r-valgrind` CI by @johnkerl in [#3318](https://github.com/single-cell-data/TileDB-SOMA/pull/3318)
- Avoid log truncation with "Last 13 lines of output" by @johnkerl in [#3313](https://github.com/single-cell-data/TileDB-SOMA/pull/3313)
- Run tests in closer to internal-dependency order by @johnkerl in [#3311](https://github.com/single-cell-data/TileDB-SOMA/pull/3311)
- Use standard name `zzz.R` for init code by @johnkerl in [#3332](https://github.com/single-cell-data/TileDB-SOMA/pull/3332)
- Address two more compiler warnings by @johnkerl in [#3320](https://github.com/single-cell-data/TileDB-SOMA/pull/3320)
- Add support for writing Seurat v5 ragged arrays
- Update docstrings for `domain` argument to `create` [#3396](https://github.com/single-cell-data/TileDB-SOMA/pull/3396)
- Vignette for new-shape feature [#3302](https://github.com/single-cell-data/TileDB-SOMA/pull/3302)
- Fix blockwise iterator + re-indexer to return re-indexed shape instead of full domain
- Docstring audit for new shape [#3399](https://github.com/single-cell-data/TileDB-SOMA/pull/3399)
- Add `check_only` support for domain/shape updates [#3400](https://github.com/single-cell-data/TileDB-SOMA/pull/3400)
- Adjust blockwise + re-indexer to return condensed matrix, not full domain [#3395](https://github.com/single-cell-data/TileDB-SOMA/pull/3395)
- \[c++\] Use core 2.27.0.rc5 [#3410](https://github.com/single-cell-data/TileDB-SOMA/pull/3410)
- \[c++\] Use core 2.27.0 [#3422](https://github.com/single-cell-data/TileDB-SOMA/pull/3422)

# tiledbsoma 1.14.5

## Changes

- Fixes a Python-only bug [#3225](https://github.com/single-cell-data/TileDB-SOMA/pull/3225)

# tiledbsoma 1.14.4

## Changes

- Add new Arrow-to-R type mapper [#3161](https://github.com/single-cell-data/TileDB-SOMA/pull/3161)
- Expose block/random writer for sparse arrays [#3204](https://github.com/single-cell-data/TileDB-SOMA/pull/3204)

# tiledbsoma 1.14.3

## Changes

- Handle `numeric` coords properly when reading arrays [3145](https://github.com/single-cell-data/TileDB-SOMA/pull/3145)

# tiledbsoma 1.14.2

## Changes

- Fixes a Python-only bug [#3074](https://github.com/single-cell-data/TileDB-SOMA/pull/3074)

# tiledbsoma 1.14.1

## Changes

- Fixes a Python-only bug [#3013](https://github.com/single-cell-data/TileDB-SOMA/pull/3013)

# tiledbsoma 1.14.0

## Changes

- New `resize` and `tiledbsoma_upgrade_shape` accessors as part of the [new-shape project](https://github.com/single-cell-data/TileDB-SOMA/issues/2407)
- Make use of timestamp ranges in libtiledbsoma
- Simplify timestamp ranges; strengthen assumptions about `tiledb_timestamp`
- Use cached timestamps in `$write()` and `$create()`
- Fix bug in blockwise iteration
- Lay groundwork for cached SOMA contexts within objects rather than re-creating contexts
- SOMA context objects are used throughout SOMA object creation
- Add value-checking for `axis` parameter when initializing blockwise reads

# tiledbsoma 1.13.1

## Changes

- Includes a fix for appended enumerations [#2903](https://github.com/single-cell-data/TileDB-SOMA/pull/2903)

## Changes

- New `resize` and `tiledbsoma_upgrade_shape` accessors as part of the [new-shape project](https://github.com/single-cell-data/TileDB-SOMA/issues/2407)
- Make use of timestamp ranges in libtiledbsoma
- Simplify timestamp ranges; strengthen assumptions about `tiledb_timestamp`
- Use cached timestamps in `$write()` and `$create()`
- Fix bug in blockwise iteration
- Lay groundwork for cached SOMA contexts within objects rather than re-creating contexts
- SOMA context objects are used throughout SOMA object creation
- Add value-checking for `axis` parameter when initializing blockwise reads

# tiledbsoma 1.13.0

## Changes

- Updates the TileDB Core dependency to 2.25.0
- The `used_shape` function is deprecated; `shape` mods are [upcoming](https://github.com/single-cell-data/TileDB-SOMA/issues/2407) scheduled for TileDB-SOMA 1.15
- Change `$reopen(mode = )` default to not flip modes; require explicit `mode` parameter to be passed
- Add `drop_levels` to `SOMAExperimentAxisQuery` -> ecosystem outgestors to drop unused factor levels

# tiledbsoma 1.12.3

## Changes

- Updates the TileDB Core dependency to 2.24.2

# tiledbsoma 1.12.2

## Changes

- No R changes; only an update for Python

# tiledbsoma 1.12.1

## Changes

- Updates to TileDB Core 2.24.1

# tiledbsoma 1.12.0

## Changes

- Updates to TileDB Core 2.24 and TileDB-R 0.28
- Connect re-indexer to blockwise iterator to return re-indexed tables and matrices

# tiledbsoma 1.11.4

## Changes

- Fixes a couple bugs observed in Python

# tiledbsoma 1.11.3

## Changes

- Fixes an intermittent segfault observed in Python

# tiledbsoma 1.11.2

## Changes

- New `reopen` method available from R

# tiledbsoma 1.11.1

## Changes

- No R changes; only an update for Python

# tiledbsoma 1.11.0

## Changes

- Add support for ingestion of `SeuratCommand` logs
- Add support for outgestion of `SeuratCommand` logs
- Add support for reading `*m` and `*p` layers from `SOMAExperimentAxisQuery`
- Add support for blockwise iteration
- Make `reopen()` a public method for all `TileDBObjects`
- Add support for resume-mode in `write_soma()`
- Push default-setting for `TileDBCreateOptions` to `$initialize()` instead of in the accessors
- Muffle warnings for missing command logs when outgesting SOMA to `Seurat`
- Have `SOMADataFrame$shape()` throw a not-yet-implemented error
- Disable running `SeuratObject::.CalcN()` when outgesting from SOMA to `Seurat`
- Clear timestamp when using `$reopen()` to reopen at the current time
- Add support for the re-indexer

# 1.10.2

- Port resume-ingest mode from Python to R

# 1.10.1

- This release contains a single Python-only bug fix

# 1.10.0

## Changes

- Add support for ingestion of `SeuratCommand` logs
- Add support for outgestion of `SeuratCommand` logs
- Add support for reading `*m` and `*p` layers from `SOMAExperimentAxisQuery`
- Add support for blockwise iteration
- Make `reopen()` a public method for all `TileDBObjects`
- Add support for resume-mode in `write_soma()`

# 1.9.5

- This release contains a single Python-only bug fix

# 1.9.4

- This release contains an in-progress blockwise iterator, work on which will be completed in a subsequent release.
- Reverted nanoarrow use to vendored versions and refactored the use throughout R and adjusted the build steps accordingly.

# 1.9.3

- This release contains a single Python-only performance improvement

# 1.9.2

- This release contains a single Python-only modification for an API backward-compatibility update

# 1.9.1

- This release contains a single Python-only modification for its build process

# 1.9.0

## Changes

- Add support for ingestion of `SeuratCommand` logs
- Add support for outgestion of `SeuratCommand` logs
- Add support for reading `*m` and `*p` layers from `SOMAExperimentAxisQuery`

# 1.8.1

## Changes

- This is a release with Python changes only

# 1.8.0

## Changes

- Add support for ingestion of `SeuratCommand` logs
- Add support for outgestion of `SeuratCommand` logs

# 1.7.0

## Changes

- Add support for registering arrays with their parent collection in `write_soma()`
- Add `write_soma.character()` to encode character vectors as SOMA arrays

# tiledbsoma 1.6.0

## Fixes

- `SOMADataFrame`, `SOMASparseNDArray`, and `SOMADenseNDArray`'s `write()` method now correctly leaves the array open in write mode

# tiledbsoma 1.5.0

## Changes

- Add support for I/O of R factors as enumerated types in `SOMADataFrame`
- Add support for writing `SummarizedExperiment` and `SingleCellExperiment` object to SOMAs
- Add support for bounding boxes for sparse arrays
- Add support for creating `SOMADataFrames` with `ordered()` columns
- Add support for I/O between SOMA and `SingleCellExperiment`
- Add support for updating `obs` and `var`
- Write string attrs as UTF-8 (Python compatibility)
- Optimize export of obsm/varm arrays to Seurat
- Add `axis_query()` method to `SOMAExperiment`
- Add `to_sparse_matrix()` method to `SOMAExperimentAxisQuery`
- Force exporting v3 assays with SeuratObject v5 installed

# tiledbsoma 1.4.0

This is the initial release of the tiledbsoma R package.

## Changes

- Added a `NEWS.md` file to track changes to the package
- `TileDBGroup` gains a `names` method to retrieve the names of group members
- Added `SOMAMeasurement` and `SOMAExperiment` classes
- spdl is now used for logging
- TileDB performance statistics can now be collected for analysis
- Added support for performing axis-based queries against a `SOMAExperiment` via the `SOMAExperimentAxisQuery` class
- `TileDBArray` class gained a `colnames()` method that returns the names of both dimensions and attributes
- Added internal helpers to centrally validate `coords` and `value_filter` arguments
- All R6 classes' `create()` method now return `self` rather than nothing
- Fixed calculating of relative paths when 1 of the URIs contains the `file://` prefix
- Added `PlatformConfig` and `SOMATileDBContext` classes to handle SOMA and TileDB configuration
- Add Seurat outgestors for `SOMAExperimentAxisQuery` objects
- Numeric coordinates passed to SOMADataFrame$read() are now automatically upcast to int64 when necessary
- Add ingestors to read data from `Seurat` objects
- Add methods for listing and accessing bundled datasets, which now includes a `SOMAExperiment` containing the pbmc_small dataset from the SeuratObject package
- New vignettes describing SOMA objects, reading data from them, and querying SOMA experiments
- Objects added to `SOMACollection`-based classes using the `add_new_*()` methods now pass through their parent context and platform config
- `SOMAExperimentAxisQuery` gained a `to_sparse_matrix()` method for retrieving data as a named sparse matrix
- `SOMAExperiment` gained `axis_query()` to construct a `SOMAExperimentAxisQuery` object
- Add SingleCellExperiment outgestor for `SOMAExperimentAxisQuery` objects
