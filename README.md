<a href="https://tiledb.com"><img src="https://github.com/TileDB-Inc/TileDB/raw/dev/doc/source/_static/tiledb-logo_color_no_margin_@4x.png" alt="TileDB logo" width="400"></a>

[![TileDB-SOMA CI](https://github.com/single-cell-data/TileDB-SOMA/actions/workflows/ci.yml/badge.svg)](https://github.com/single-cell-data/TileDB-SOMA/actions/workflows/ci.yml)
[![TileDB-SOMA R CI](https://github.com/single-cell-data/TileDB-SOMA/actions/workflows/r-ci.yml/badge.svg)](https://github.com/single-cell-data/TileDB-SOMA/actions/workflows/r-ci.yml)
[![codecov](https://codecov.io/github/single-cell-data/TileDB-SOMA/branch/main/graph/badge.svg)](https://codecov.io/github/single-cell-data/TileDB-SOMA)


# TileDB-SOMA

[SOMA](https://github.com/single-cell-data/SOMA/tree/main) – for “Stack Of Matrices, Annotated” – is a flexible, extensible, and open-source API enabling access to data in a variety of formats. The driving use case of SOMA is for single-cell data in the form of annotated matrices where observations are frequently cells and features are genes, proteins, or genomic regions.

The TileDB-SOMA package is a C++ library with APIs in Python and R, using [TileDB
Embedded](https://github.com/TileDB-Inc/TileDB) to implement the
[SOMA specification](https://github.com/single-cell-data/SOMA/blob/main/abstract_specification.md).

Get started on using TileDB-SOMA:

* [Quick start](#quick-start).
* Python [documentation](#tiledb-soma) | [tutorials](#tiledb-soma). *Under development*.
* R [documentation](#tiledb-soma) | [tutorials](#tiledb-soma). *Under development*.

## What can TileDB-SOMA do?

Intended to be used for single-cell data, TileDB-SOMA provides Python and R APIs to allow for storage and data access patterns at scale and for larger-than-memory operations:

* Create and write large volumes of data.
* Open and read data at low latency, locally and from the cloud.
* Query and access interconnected arrays efficiently and at low latency.

TiledDB-SOMA provides interoperability with existing single-cell toolkits:

* Load and create [anndata](https://anndata.readthedocs.io/en/latest/) objects.
* Load and create [Seurat](https://anndata.readthedocs.io/en/latest/) objects. *Coming soon*.

TiledDB-SOMA provides interoperability with existing Python or R data structures:

* From Python create PyArrow objects, SciPy sparse matrices, NumPy arrays, and pandas data frames.
* From R create R Arrow objects, sparse matrices (via the [Matrix](https://cran.r-project.org/package=Matrix) package), and standard data frames and (dense) matrices.


## Community

* Please join the [TileDB Slack community](https://tiledb-community.slack.com/join/shared_invite/zt-ndq1ipwl-QcithaWG6j1BImtuQGSpag#/shared-invite/email) with dedicated channel `#genomics`
* Please join the [CZI Slack community](https://cziscience.slack.com/), with dedicated
channel `#cell-census-users`


## Quick start

### Documentation

The TileDB-SOMA doc-site (under development), contains the reference documentation and tutorials.

Reference documentation can also be accessed directly from Python `help(tiledsoma)` or R `help(package = "tiledbsoma")`.

### Main SOMA Objects

The capabilities of TileDB-SOMA lay on the different read, access, and query patterns that each of the main implementations of SOMA objects provide:

* `DenseNDArray` a dense, N-dimensional array, with offset (zero-based) integer indexing on each dimension.
* `SparseNDArray` same as `DenseNDArray` but sparse, and it supports point indexing (disjoint index access)
* `DataFrame`  a multi-column table with a user-defined columns names and value types, with support for point indexing. 
* `Collection` a persistent container of named SOMA objects.
* `Experiment` is a class that represents a single-cell experiment. It always contains two objects:
	* `obs`: a  `DataFrame` with primary annotations on the observation axis.
	* `ms`: a  `Collection` of measurements, each composed of `X` matrices and axis annotation matrices or data frames.

### APIs quick start

* [Python quick start](https://github.com/single-cell-data/TileDB-SOMA/wiki/Python-quick-start)
* [R quick start](https://github.com/single-cell-data/TileDB-SOMA/wiki/R-quick-start)

## Issues and contacts

* Any/all questions, comments, and concerns are welcome at the [GitHub new-issue page](https://github.com/single-cell-data/TileDB-SOMA/issues/new/choose) -- or, you can also browse [existing issues](https://github.com/single-cell-data/TileDB-SOMA/issues)
* If you believe you have found a security issue, in lieu of filing an issue please responsibly disclose it by contacting [security@tiledb.com](mailto:security@tiledb.com)

## Branches

This branch, `main`, implements the [updated specfication](https://github.com/single-cell-data/SOMA/blob/main/abstract_specification.md).  Please also see the `main-old` branch which implements the [original specification](https://github.com/single-cell-data/TileDB-SOMA/blob/main-old/spec/specification.md).

## Developer information

* [TileDB-SOMA wiki](https://github.com/single-cell-data/TileDB-SOMA/wiki)
* [Build instructions for developers](libtiledbsoma/README.md)

## Code of Conduct

All participants in TileDB spaces are expected to adhere to high standards of
professionalism in all interactions. This repository is governed by the
specific standards and reporting procedures detailed in depth in the
[TileDB core repository Code Of Conduct](
https://github.com/TileDB-Inc/TileDB/blob/dev/CODE_OF_CONDUCT.md).

<!-- links -->
[tiledb]: https://github.com/TileDB-Inc/TileDB
