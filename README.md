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
* Python [documentation](.) | [tutorials](.). *Under development*.
* R [documentation](.) | [tutorials](.). *Under development*.

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
* From R create RArrow objects, Matrix sparse matrices, and R-base dense matrices.


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

### Python quick start

#### Installation

TileDB-SOMA is available on [PyPI](https://pypi.org/project/tiledbsoma/) and can be installed via `pip` as indicated below. Full installation instructions can be found [here](https://github.com/single-cell-data/TileDB-SOMA/blob/main/apis/python/README.md).

```bash
pip install tiledbsoma
```

#### Usage examples

##### Building a SOMA object

SOMA objects can be created with their respective `create()` methods and then need to be populated in specific ways depending on their types. 

However, a [`SOMAExperiment`](https://github.com/single-cell-data/SOMA/blob/main/abstract_specification.md#somaexperiment-and-somameasurement) can be easily created from and anndata object or a `*h5ad` file. Here, one is created from a `*.h5ad` file.

```python
import tiledbsoma.io

# Create a and write a SOMA Experiment, source data https://github.com/chanzuckerberg/cellxgene/raw/main/example-dataset/pbmc3k.h5ad
pbmc3k_uri = tiledbsoma.io.from_h5ad("./pbmc3k", input_path = "pbmc3k.h5ad", measurement_name = "RNA")
``` 

##### Reading and querying SOMA objects

SOMA objects can be opened using `tildedbsoma.open()`. 

The contents of `DataFrame`, `SparseNDArray` and `DenseNDArray` can be accessed with their respective `read()` methods. For  `DataFrame` and `SparseNDArray` the method returns an iterator useful for larger-than-memory operations.

For example you can open the `SOMAExperiment` created above and then read the contents of `obs` which is a `SOMADataFrame`. 

In addition, this example shows how you can query for observations with `louvian` values of 'Megakaryocytes' and 'CD4 T cells', and `n_genes` greater than 500.

```python
import tiledbsoma

with tiledbsoma.open(pbmc3k_uri) as pbmc3k_soma:
    pbmc3k_obs_slice = pbmc3k_soma.obs.read(
        value_filter="n_genes >500 and louvain in ['Megakaryocytes', 'CD4 T cells']"
    )
    
    # Concatenate iterator to pyarrow.Table
    pbmc3k_obs_slice.concat()
```

The result is a `pyarrow.Table` containing a slice based on the specified filters.

```bash
pyarrow.Table
soma_joinid: int64
obs_id: large_string
n_genes: int64
percent_mito: float
n_counts: float
louvain: large_string
----
soma_joinid: [[0,2,8,11,12,...,2617,2621,2626,2631,2637]]
obs_id: [["AAACATACAACCAC-1","AAACATTGATCAGC-1","AAACGCTGTAGCCA-1","AAACTTGATCCAGA-1","AAAGAGACGAGATA-1",...,"TTGTAGCTAGCTCA-1","TTTAGCTGATACCG-1","TTTCACGAGGTTCA-1","TTTCCAGAGGTGAG-1","TTTGCATGCCTCAC-1"]]
n_genes: [[781,1131,533,751,866,...,933,887,721,873,724]]
percent_mito: [[0.030177759,0.008897362,0.011764706,0.010887772,0.010788382,...,0.02224871,0.022875817,0.013261297,0.0068587107,0.008064516]]
n_counts: [[2419,3147,1275,2388,2410,...,2517,2754,2036,2187,1984]]
louvain: [["CD4 T cells","CD4 T cells","CD4 T cells","CD4 T cells","CD4 T cells",...,"CD4 T cells","CD4 T cells","CD4 T cells","CD4 T cells","CD4 T cells"]]
```

##### Iterators for larger-than-memory operations

As stated above the `read()` methods of `DataFrame` and `SparseNDArray` return an iterator. The batch size can be specified a in the `soma.init_buffer_bytes` config option, for this is example it is set to 100 Bytes:

```python
context = tiledbsoma.options.SOMATileDBContext()
context = context.replace(tiledb_config = {"soma.init_buffer_bytes": 100})

with tiledbsoma.open(pbmc3k_uri, context = context) as pbmc3k_soma:
    
    pbmc3k_obs = pbmc3k_soma.obs.read()
  
    counter = 1
    for pbmc3k_obs_chunk in pbmc3k_obs:
        
        # Perform operations
        # pbmc3k_obs_chunk is a pyArrow.Table
        
        counter += 1

print(counter)
```

The counter indicates the number of iterations performed

```bash
441
```

### R quick start

#### Installation

TileDB-SOMA is available on  R-universe and can be installed directly from R as indicated below. Full installation instructions can be found [here](https://github.com/single-cell-data/TileDB-SOMA/blob/main/apis/r/README.md).

```r
install.packages('tiledbsoma', repos = c('https://tiledb-inc.r-universe.dev', 'https://cloud.r-project.org'))
```

#### Usage examples

*Coming Soon*

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
