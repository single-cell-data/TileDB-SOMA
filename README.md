<a href="https://tiledb.com"><img src="https://github.com/TileDB-Inc/TileDB/raw/dev/doc/source/_static/tiledb-logo_color_no_margin_@4x.png" alt="TileDB logo" width="400"></a>

[![TileDB-SOMA CI](https://github.com/single-cell-data/TileDB-SOMA/actions/workflows/ci.yml/badge.svg)](https://github.com/single-cell-data/TileDB-SOMA/actions/workflows/ci.yml)
[![TileDB-SOMA R CI](https://github.com/single-cell-data/TileDB-SOMA/actions/workflows/r-ci.yml/badge.svg)](https://github.com/single-cell-data/TileDB-SOMA/actions/workflows/r-ci.yml)
[![PyPI version](https://badge.fury.io/py/tiledbsoma.svg)](https://badge.fury.io/py/tiledbsoma)
[![tiledbsoma status badge](https://tiledb-inc.r-universe.dev/badges/tiledbsoma)](https://tiledb-inc.r-universe.dev)
[![codecov](https://codecov.io/github/single-cell-data/TileDB-SOMA/branch/main/graph/badge.svg)](https://codecov.io/github/single-cell-data/TileDB-SOMA)


# TileDB-SOMA

[SOMA](https://github.com/single-cell-data/SOMA/tree/main) – for “Stack Of Matrices, Annotated” – is a flexible, extensible, and open-source API enabling access to data in a variety of formats. The driving use case of SOMA is for single-cell data in the form of annotated matrices where observations are frequently cells and features are genes, proteins, or genomic regions.

The TileDB-SOMA package is a C++ library with APIs in Python and R, using [TileDB
Embedded](https://github.com/TileDB-Inc/TileDB) to implement the
[SOMA specification](https://github.com/single-cell-data/SOMA/blob/main/abstract_specification.md).

Get started on using TileDB-SOMA:

* [Quick start](#quick-start).
* Python [documentation](https://tiledbsoma.readthedocs.io/en/latest/python-api.html). 
* R [documentation](https://single-cell-data.github.io/TileDB-SOMA/).


## What Can TileDB-SOMA Do?

Intended to be used for single-cell data, TileDB-SOMA provides Python and R APIs to allow for storage and data access patterns at scale and for larger-than-memory operations:


* Create and write large volumes of data.
* Open and read data at low latency, locally and from the cloud.
* Query and access interconnected arrays efficiently and at low latency.

TileDB-SOMA provides interoperability with existing single-cell toolkits:

* Load and create [AnnData](https://anndata.readthedocs.io/en/latest/) objects.
* Load and create [Seurat](https://satijalab.org/seurat/) objects. *Coming soon*.

TileDB-SOMA provides interoperability with existing Python or R data structures:

* From Python create PyArrow objects, SciPy sparse matrices, NumPy arrays, and pandas data frames.
* From R create R Arrow objects, sparse matrices (via the [Matrix](https://cran.r-project.org/package=Matrix) package), and standard data frames and (dense) matrices.


## Community

* Please join the [TileDB Slack community](https://tiledb-community.slack.com/join/shared_invite/zt-ndq1ipwl-QcithaWG6j1BImtuQGSpag#/shared-invite/email) with dedicated channel `#genomics`.
* Please join the [CZI Slack community](https://cziscience.slack.com/join/shared_invite/zt-czl1kp2v-sgGpY4RxO3bPYmFg2XlbZA#/shared-invite/email), with dedicated
channel `#cellxgene-census-users`.


## APIs Installation and Quick Start

* [Python installation and quick start](https://github.com/single-cell-data/TileDB-SOMA/wiki/Python-quick-start)
* [R installation and quick start](https://github.com/single-cell-data/TileDB-SOMA/wiki/R-quick-start)

## API Documentation

The TileDB-SOMA doc-site ([Python](https://tiledbsoma.readthedocs.io/en/latest/python-api.html)|[R](https://single-cell-data.github.io/TileDB-SOMA/)), contains the reference documentation and tutorials.

Reference documentation can also be accessed directly from Python `help(tiledsoma)` or R `help(package = "tiledbsoma")`.

## Main SOMA Objects

The capabilities of TileDB-SOMA lay on the different read, access, and query patterns that each of the main implementations of SOMA objects provide:

* `DenseNDArray` is a dense, N-dimensional array, with offset (zero-based) integer indexing on each dimension.
* `SparseNDArray` is the same as `DenseNDArray` but sparse, and supports point indexing (disjoint index access).
* `DataFrame` is a multi-column table with a user-defined columns names and value types, with support for point indexing. 
* `Collection` is a persistent container of named SOMA objects.
* `Experiment` is a class that represents a single-cell experiment. It always contains two objects:
	* `obs`: a  `DataFrame` with primary annotations on the observation axis.
	* `ms`: a  `Collection` of measurements, each composed of `X` matrices and axis annotation matrices or data frames (e.g. `var`, `varm`, `obsm`, etc).

## Who Is Using SOMA?

* [CZ CELLxGENE Discover](https://cellxgene.cziscience.com/) to build its [Census](https://github.com/chanzuckerberg/cellxgene-census/), which provides efficient access and querying to a corpus containing nearly 50 million cells, compiled from 700+ datasets.

If you are interested in listing any projects here please contact us at [soma@chanzuckerberg.com](mailto:soma@chanzuckerberg.com).


## Issues and Contacts

* Any/all questions, comments, and concerns are welcome at the [GitHub new-issue page](https://github.com/single-cell-data/TileDB-SOMA/issues/new/choose) -- or, you can also browse [existing issues](https://github.com/single-cell-data/TileDB-SOMA/issues).
* If you believe you have found a security issue, in lieu of filing an issue please responsibly disclose it by contacting [security@tiledb.com](mailto:security@tiledb.com).


## Branches

This branch, `main`, implements the [updated specfication](https://github.com/single-cell-data/SOMA/blob/main/abstract_specification.md).  Please also see the `main-old` branch which implements the [original specification](https://github.com/single-cell-data/TileDB-SOMA/blob/main-old/spec/specification.md).


## Developer Information

* [TileDB-SOMA wiki](https://github.com/single-cell-data/TileDB-SOMA/wiki).
* [Build instructions for developers](libtiledbsoma/README.md).


## Code of Conduct

All participants in TileDB spaces are expected to adhere to high standards of
professionalism in all interactions. This repository is governed by the
specific standards and reporting procedures detailed in depth in the
[TileDB core repository Code Of Conduct](
https://github.com/TileDB-Inc/TileDB/blob/dev/CODE_OF_CONDUCT.md).

<!-- links -->
[tiledb]: https://github.com/TileDB-Inc/TileDB
