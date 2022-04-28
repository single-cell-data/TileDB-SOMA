This is test code for reading ANN data and writing into a TileDB nested group structure.

# TL;DR

```
./desc-ann.py ./anndata/pbmc3k_processed.h5ad

./ingestor.py ./anndata/pbmc3k_processed.h5ad ./tiledb-data/pbmc3k_processed

./desc-tiledb.py ./tiledb-data/pbmc3k_processed
```

# Installation

This requires [`tiledb`](https://github.com/TileDB-Inc/TileDB-Py) 0.14.1 or above, in addition to other dependencies
in [setup.cfg](./setup.cfg).

After `cd` to `apis/python`:

```
pip install --editable .
```

Optionally, if you prefer, instead:

```
python -m venv venv
. ./venv/bin/activate
pip install .
```

Then:

```
python -m pytest tests
```

# Overview

* Sample data:
  * [anndata](./anndata) contains some files from [https://cellxgene.cziscience.com](https://cellxgene.cziscience.com)
  * The most important reference is `anndata/pbmc3k_processed.h5ad`
* Code:
  * [./src/tiledbsc](./src/tiledbsc)
* Inspecting HDF5 input files
  * `./desc-ann.py ./anndata/pbmc3k_processed.h5ad`
* Ingesting
  * `./ingestor.py ./anndata/pbmc3k_processed.h5ad`
  * Output is in `tiledb-data/pbmc3k_processed`
  * Cloud-upload test:
    * `ingestor.py ./anndata/pbmc3k_processed.h5ad tiledb://johnkerl-tiledb/s3://tiledb-johnkerl/wpv2-test-001`
* Inspecting TileDB output groups
  * `./desc-tiledb.py ./tiledb-data/pbmc3k_processed`

# Status

* Most important: input files beyond `./anndata/pbmc3k_processed.h5ad` need to be validated
  * Processing of this datafile is validated in this README, unit-tests for this package, and the CI job for this repo
* The `uns` arrays from HDF5 files are currently not processed
* Dimensioning of `X/data` and `X/raw` is non-homogeneous in cellxgene datasets looked at so far (in contrast to Seurat data) -- this needs discussion with CZI
* Labeling in `obs`, `var`, etc needs to be more aligned with that used by the [tiledbsc R package](https://github.com/TileDB-Inc/tiledbsc).

# Details

## Expected input format

`./desc-ann.py ./anndata/pbmc3k_processed.h5ad`

<details>

```
================================================================ ./anndata/pbmc3k_processed.h5ad
ANNDATA SUMMARY:
AnnData object with n_obs × n_vars = 2638 × 1838
    obs: 'n_genes', 'percent_mito', 'n_counts', 'louvain'
    var: 'n_cells'
    uns: 'draw_graph', 'louvain', 'louvain_colors', 'neighbors', 'pca', 'rank_genes_groups'
    obsm: 'X_pca', 'X_tsne', 'X_umap', 'X_draw_graph_fr'
    varm: 'PCs'
    obsp: 'distances', 'connectivities'
X IS A    <class 'numpy.ndarray'>
  X SHAPE   (2638, 1838)
  OBS  LEN  2638
  VAR  LEN  1838
RAW X IS A    <class 'scipy.sparse._csr.csr_matrix'>
  X SHAPE   (2638, 13714)
OBS IS A <class 'pandas.core.frame.DataFrame'>
  OBS  KEYS ['n_genes', 'percent_mito', 'n_counts', 'louvain']
VAR IS A <class 'pandas.core.frame.DataFrame'>
  VAR  KEYS ['n_cells']
OBSM KEYS ['X_pca', 'X_tsne', 'X_umap', 'X_draw_graph_fr']
  OBSM X_pca IS A <class 'numpy.ndarray'>
  OBSM X_tsne IS A <class 'numpy.ndarray'>
  OBSM X_umap IS A <class 'numpy.ndarray'>
  OBSM X_draw_graph_fr IS A <class 'numpy.ndarray'>
VARM KEYS ['PCs']
  VARM PCs IS A <class 'numpy.ndarray'>
OBSP KEYS ['distances', 'connectivities']
  OBSP distances IS A <class 'scipy.sparse.csr.csr_matrix'>
  OBSP connectivities IS A <class 'scipy.sparse.csr.csr_matrix'>
VARP KEYS []
```

See also:

```
h5ls     anndata/pbmc3k_processed.h5ad
h5ls -r  anndata/pbmc3k_processed.h5ad
h5ls -vr anndata/pbmc3k_processed.h5ad
```

</details>

## Expected ingestion progress

`./ingestor.py ./anndata/pbmc3k_processed.h5ad ./tiledb-data/pbmc3k_processed`

<details>

```
START  SOMA.from_h5ad ./anndata/pbmc3k_processed.h5ad -> ./tiledb-data/pbmc3k_processed
  START  READING ./anndata/pbmc3k_processed.h5ad
  FINISH READING ./anndata/pbmc3k_processed.h5ad
  START  DECATEGORICALIZING
  FINISH DECATEGORICALIZING
  START  WRITING ./tiledb-data/pbmc3k_processed
    START  WRITING ./tiledb-data/pbmc3k_processed/X/data
    FINISH WRITING ./tiledb-data/pbmc3k_processed/X/data
    START  WRITING ./tiledb-data/pbmc3k_processed/X/raw
    FINISH WRITING ./tiledb-data/pbmc3k_processed/X/raw
    START  WRITING ./tiledb-data/pbmc3k_processed/obs
    FINISH WRITING ./tiledb-data/pbmc3k_processed/obs
    START  WRITING ./tiledb-data/pbmc3k_processed/var
    FINISH WRITING ./tiledb-data/pbmc3k_processed/var
    START  WRITING ./tiledb-data/pbmc3k_processed/obsm/X_pca
    FINISH WRITING ./tiledb-data/pbmc3k_processed/obsm/X_pca
    START  WRITING ./tiledb-data/pbmc3k_processed/obsm/X_tsne
    FINISH WRITING ./tiledb-data/pbmc3k_processed/obsm/X_tsne
    START  WRITING ./tiledb-data/pbmc3k_processed/obsm/X_umap
    FINISH WRITING ./tiledb-data/pbmc3k_processed/obsm/X_umap
    START  WRITING ./tiledb-data/pbmc3k_processed/obsm/X_draw_graph_fr
    FINISH WRITING ./tiledb-data/pbmc3k_processed/obsm/X_draw_graph_fr
    START  WRITING ./tiledb-data/pbmc3k_processed/varm/PCs
    FINISH WRITING ./tiledb-data/pbmc3k_processed/varm/PCs
    START  WRITING ./tiledb-data/pbmc3k_processed/obsp/distances
    FINISH WRITING ./tiledb-data/pbmc3k_processed/obsp/distances
    START  WRITING ./tiledb-data/pbmc3k_processed/obsp/connectivities
    FINISH WRITING ./tiledb-data/pbmc3k_processed/obsp/connectivities
  FINISH WRITING ./tiledb-data/pbmc3k_processed
FINISH SOMA.from_h5ad ./anndata/pbmc3k_processed.h5ad -> ./tiledb-data/pbmc3k_processed
```

</details>

## Expected output format

`./desc-tiledb.py ./tiledb-data/pbmc3k_processed`

<details>

```
================================================================
X/data:
keys ['data', 'obs_id', 'var_id']
OrderedDict([('data', array([-0.13904382, -0.0708308 , -0.54010755, ..., -0.30667225,
       -0.15601887,  3.3442626 ])), ('obs_id', array([b'AAACATACAACCAC-1', b'AAACATACAACCAC-1', b'AAACATACAACCAC-1', ...,
       b'TTTGCATGCCTCAC-1', b'TTTGCATGCCTCAC-1', b'TTTGCATGCCTCAC-1'],
      dtype=object)), ('var_id', array([b'AAGAB', b'AAR2', b'AATF', ..., b'ZUFSP', b'ZWINT', b'ZYX'],
      dtype=object))])
ArraySchema(
  domain=Domain(*[
    Dim(name='obs_id', domain=(None, None), tile=None, dtype='|S0', var=True),
    Dim(name='var_id', domain=(None, None), tile=None, dtype='|S0', var=True),
  ]),
  attrs=[
    Attr(name='data', dtype='float64', var=False, nullable=False),
  ],
  cell_order='row-major',
  tile_order='row-major',
  capacity=10000,
  sparse=True,
  allows_duplicates=False,
)

X/raw:
keys ['raw', 'obs_id', 'var_id']
OrderedDict([('raw', array([0.        , 0.        , 1.60943794, ..., 0.        , 0.        ,
       0.        ])), ('obs_id', array([b'AAACATACAACCAC-1', b'AAACATTGAGCTAC-1', b'AAACATTGATCAGC-1', ...,
       b'TTTGCATGAGAGGC-1', b'TTTGCATGCCTCAC-1', b'TTTGCATGCCTCAC-1'],
      dtype=object)), ('var_id', array([b'7SK-2', b'7SK-2', b'7SK-2', ..., b'hsa-mir-8072',
       b'hsa-mir-1199', b'hsa-mir-8072'], dtype=object))])
ArraySchema(
  domain=Domain(*[
    Dim(name='obs_id', domain=(None, None), tile=None, dtype='|S0', var=True),
    Dim(name='var_id', domain=(None, None), tile=None, dtype='|S0', var=True),
  ]),
  attrs=[
    Attr(name='raw', dtype='float64', var=False, nullable=False),
  ],
  cell_order='row-major',
  tile_order='row-major',
  capacity=10000,
  sparse=True,
  allows_duplicates=False,
)

----------------------------------------------------------------
obs:
keys ['n_genes', 'percent_mito', 'n_counts', 'louvain', 'index']
ArraySchema(
  domain=Domain(*[
    Dim(name='index', domain=(None, None), tile=None, dtype='|S0', var=True, filters=FilterList([ZstdFilter(level=-1), ])),
  ]),
  attrs=[
    Attr(name='n_genes', dtype='int64', var=False, nullable=False, filters=FilterList([ZstdFilter(level=-1), ])),
    Attr(name='percent_mito', dtype='float32', var=False, nullable=False, filters=FilterList([ZstdFilter(level=-1), ])),
    Attr(name='n_counts', dtype='float32', var=False, nullable=False, filters=FilterList([ZstdFilter(level=-1), ])),
    Attr(name='louvain', dtype='<U0', var=True, nullable=False, filters=FilterList([ZstdFilter(level=-1), ])),
  ],
  cell_order='row-major',
  tile_order='row-major',
  capacity=10000,
  sparse=True,
  allows_duplicates=False,
)

----------------------------------------------------------------
var:
keys ['n_cells', 'index']
ArraySchema(
  domain=Domain(*[
    Dim(name='index', domain=(None, None), tile=None, dtype='|S0', var=True, filters=FilterList([ZstdFilter(level=-1), ])),
  ]),
  attrs=[
    Attr(name='n_cells', dtype='int64', var=False, nullable=False, filters=FilterList([ZstdFilter(level=-1), ])),
  ],
  cell_order='row-major',
  tile_order='row-major',
  capacity=10000,
  sparse=True,
  allows_duplicates=False,
)


----------------------------------------------------------------
obsm:
   file:///Users/johnkerl/git/johnkerl/TileDB-SingleCell/util/tiledb-data/pbmc3k_processed/obsm/X_umap (2638, 2)
ArraySchema(
  domain=Domain(*[
    Dim(name='__dim_0', domain=(0, 2637), tile=2638, dtype='uint64'),
    Dim(name='__dim_1', domain=(0, 1), tile=2, dtype='uint64'),
  ]),
  attrs=[
    Attr(name='', dtype='float64', var=False, nullable=False),
  ],
  cell_order='row-major',
  tile_order='row-major',
  capacity=10000,
  sparse=False,
)

   file:///Users/johnkerl/git/johnkerl/TileDB-SingleCell/util/tiledb-data/pbmc3k_processed/obsm/X_tsne (2638, 2)
ArraySchema(
  domain=Domain(*[
    Dim(name='__dim_0', domain=(0, 2637), tile=2638, dtype='uint64'),
    Dim(name='__dim_1', domain=(0, 1), tile=2, dtype='uint64'),
  ]),
  attrs=[
    Attr(name='', dtype='float64', var=False, nullable=False),
  ],
  cell_order='row-major',
  tile_order='row-major',
  capacity=10000,
  sparse=False,
)

   file:///Users/johnkerl/git/johnkerl/TileDB-SingleCell/util/tiledb-data/pbmc3k_processed/obsm/X_draw_graph_fr (2638, 2)
ArraySchema(
  domain=Domain(*[
    Dim(name='__dim_0', domain=(0, 2637), tile=2638, dtype='uint64'),
    Dim(name='__dim_1', domain=(0, 1), tile=2, dtype='uint64'),
  ]),
  attrs=[
    Attr(name='', dtype='float64', var=False, nullable=False),
  ],
  cell_order='row-major',
  tile_order='row-major',
  capacity=10000,
  sparse=False,
)

   file:///Users/johnkerl/git/johnkerl/TileDB-SingleCell/util/tiledb-data/pbmc3k_processed/obsm/X_pca (2638, 50)
ArraySchema(
  domain=Domain(*[
    Dim(name='__dim_0', domain=(0, 2637), tile=2638, dtype='uint64'),
    Dim(name='__dim_1', domain=(0, 49), tile=50, dtype='uint64'),
  ]),
  attrs=[
    Attr(name='', dtype='float32', var=False, nullable=False),
  ],
  cell_order='row-major',
  tile_order='row-major',
  capacity=10000,
  sparse=False,
)


----------------------------------------------------------------
varm:
   file:///Users/johnkerl/git/johnkerl/TileDB-SingleCell/util/tiledb-data/pbmc3k_processed/varm/PCs (1838, 50)
ArraySchema(
  domain=Domain(*[
    Dim(name='__dim_0', domain=(0, 1837), tile=1838, dtype='uint64'),
    Dim(name='__dim_1', domain=(0, 49), tile=50, dtype='uint64'),
  ]),
  attrs=[
    Attr(name='', dtype='float32', var=False, nullable=False),
  ],
  cell_order='row-major',
  tile_order='row-major',
  capacity=10000,
  sparse=False,
)


----------------------------------------------------------------
obsp:
   file:///Users/johnkerl/git/johnkerl/TileDB-SingleCell/util/tiledb-data/pbmc3k_processed/obsp/connectivities (2638, 2638)
ArraySchema(
  domain=Domain(*[
    Dim(name='__dim_0', domain=(0, 2637), tile=2638, dtype='uint64'),
    Dim(name='__dim_1', domain=(0, 2637), tile=2638, dtype='uint64'),
  ]),
  attrs=[
    Attr(name='', dtype='float64', var=False, nullable=False),
  ],
  cell_order='row-major',
  tile_order='row-major',
  capacity=10000,
  sparse=False,
)

   file:///Users/johnkerl/git/johnkerl/TileDB-SingleCell/util/tiledb-data/pbmc3k_processed/obsp/distances (2638, 2638)
ArraySchema(
  domain=Domain(*[
    Dim(name='__dim_0', domain=(0, 2637), tile=2638, dtype='uint64'),
    Dim(name='__dim_1', domain=(0, 2637), tile=2638, dtype='uint64'),
  ]),
  attrs=[
    Attr(name='', dtype='float64', var=False, nullable=False),
  ],
  cell_order='row-major',
  tile_order='row-major',
  capacity=10000,
  sparse=False,
)

```

</details>

See also:

```
import tiledb
print(tiledb.group.Group('tiledb-data/pbmc3k_processed')._dump(True))
```

## Diversity of formats in HDF5 files

Due to file-size restrictions on GitHub, not all the following are cached in this repo.
Nonetheless this is a selection from
[https://cellxgene.cziscience.com](https://cellxgene.cziscience.com), as well as some
raw-sensor data (`subset_100_100`):

<details>

```
for x in $(ls -Sr ~/ann/*.h5ad); do
  echo ================================================== $x
   ./desc-ann.py $x | grep IS.A
done
```

```
================================================== /Users/johnkerl/ann/subset_100_100.h5ad
X IS A    <class 'scipy.sparse._csr.csr_matrix'>
OBS IS A <class 'pandas.core.frame.DataFrame'>
VAR IS A <class 'pandas.core.frame.DataFrame'>
================================================== /Users/johnkerl/ann/pbmc3k_processed.h5ad
X IS A    <class 'numpy.ndarray'>
OBS IS A <class 'pandas.core.frame.DataFrame'>
VAR IS A <class 'pandas.core.frame.DataFrame'>
  OBSM X_pca IS A <class 'numpy.ndarray'>
  OBSM X_tsne IS A <class 'numpy.ndarray'>
  OBSM X_umap IS A <class 'numpy.ndarray'>
  OBSM X_draw_graph_fr IS A <class 'numpy.ndarray'>
  VARM PCs IS A <class 'numpy.ndarray'>
  OBSP distances IS A <class 'scipy.sparse._csr.csr_matrix'>
  OBSP connectivities IS A <class 'scipy.sparse._csr.csr_matrix'>
================================================== /Users/johnkerl/ann/local3.h5ad
X IS A    <class 'scipy.sparse._csr.csr_matrix'>
OBS IS A <class 'pandas.core.frame.DataFrame'>
VAR IS A <class 'pandas.core.frame.DataFrame'>
  OBSM X_umap IS A <class 'numpy.ndarray'>
================================================== /Users/johnkerl/ann/human-kidney-tumors-wilms.h5ad
X IS A    <class 'scipy.sparse._csr.csr_matrix'>
OBS IS A <class 'pandas.core.frame.DataFrame'>
VAR IS A <class 'pandas.core.frame.DataFrame'>
  OBSM X_pca IS A <class 'numpy.ndarray'>
  OBSM X_umap IS A <class 'numpy.ndarray'>
================================================== /Users/johnkerl/ann/longitudinal-profiling-49.h5ad
X IS A    <class 'scipy.sparse._csr.csr_matrix'>
OBS IS A <class 'pandas.core.frame.DataFrame'>
VAR IS A <class 'pandas.core.frame.DataFrame'>
  OBSM X_umap IS A <class 'numpy.ndarray'>
================================================== /Users/johnkerl/ann/single-cell-transcriptomes.h5ad
X IS A    <class 'scipy.sparse._csr.csr_matrix'>
OBS IS A <class 'pandas.core.frame.DataFrame'>
VAR IS A <class 'pandas.core.frame.DataFrame'>
  OBSM X_pca IS A <class 'numpy.ndarray'>
  OBSM X_umap IS A <class 'numpy.ndarray'>
================================================== /Users/johnkerl/ann/vieira19_Alveoli_and_parenchyma_anonymised.processed.h5ad
X IS A    <class 'scipy.sparse._csc.csc_matrix'>
OBS IS A <class 'pandas.core.frame.DataFrame'>
VAR IS A <class 'pandas.core.frame.DataFrame'>
  OBSM X_umap_hm IS A <class 'numpy.ndarray'>
  VARM PCs IS A <class 'numpy.ndarray'>
================================================== /Users/johnkerl/ann/local2.h5ad
X IS A    <class 'scipy.sparse._csr.csr_matrix'>
OBS IS A <class 'pandas.core.frame.DataFrame'>
VAR IS A <class 'pandas.core.frame.DataFrame'>
  OBSM X_umap IS A <class 'numpy.ndarray'>
================================================== /Users/johnkerl/ann/acute-covid19-cohort.h5ad
X IS A    <class 'scipy.sparse._csr.csr_matrix'>
OBS IS A <class 'pandas.core.frame.DataFrame'>
VAR IS A <class 'pandas.core.frame.DataFrame'>
  OBSM X_umap IS A <class 'numpy.ndarray'>
================================================== /Users/johnkerl/ann/autoimmunity-pbmcs.h5ad
X IS A    <class 'scipy.sparse._csr.csr_matrix'>
OBS IS A <class 'pandas.core.frame.DataFrame'>
VAR IS A <class 'pandas.core.frame.DataFrame'>
  OBSM X_umap IS A <class 'numpy.ndarray'>
```

See also:

```
for x in ~/ann/*.h5ad; do
  echo ================================================ $x
  h5ls $x
done
```

</details>

# Notes

* `os.path.join` is used here but may not be appropriate if this package is run on Windows:
  * `/` has been accepted in Windows paths for some years now
  * `\` is not accepted for forming URIs
  * So, perhaps safer would be to always join on `/` regardless of platform.
