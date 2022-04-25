This is test code for reading ANN data and writing into a TileDB nested group structure.

# TL;DR

```
./desc-ann.py ./anndata/pbmc3k_processed.h5ad

./ingestor.py ./anndata/pbmc3k_processed.h5ad ./tiledb-data/pbmc3k_processed

./desc-tiledb.py ./tiledb-data/pbmc3k_processed
```

# Overview

* Sample data:
  * [anndata](./anndata) contains some files from [https://cellxgene.cziscience.com](https://cellxgene.cziscience.com)
  * The most important reference is `anndata/pbmc3k_processed.h5ad`
* Code:
  * [mSCGroup.py](./mSCGroup.py)
* Inspecting HD5 input files
  * `./desc-ann.py ./anndata/pbmc3k_processed.h5ad`
* Ingesting
  * `./ingestor.py ./anndata/pbmc3k_processed.h5ad`
  * Output is in `tiledb-data/pbmc3k_processed`
  * Cloud-upload test:
    * `ingestor.py ./anndata/pbmc3k_processed.h5ad tiledb://johnkerl-tiledb/s3://tiledb-johnkerl/wpv2-test-001`
* Inspecting TileDB output groups
  * `./desc-tiledb.py ./tiledb-data/pbmc3k_processed`

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

</details>

## Expected ingestion progress

`./ingestor.py ./anndata/pbmc3k_processed.h5ad ./tiledb-data/pbmc3k_processed`

<details>

```
START  SCGroup.from_h5ad ./anndata/pbmc3k_processed.h5ad -> ./tiledb-data/pbmc3k_processed
  START  READING ./anndata/pbmc3k_processed.h5ad
  FINISH READING ./anndata/pbmc3k_processed.h5ad
  START  DECATEGORICALIZING
  FINISH DECATEGORICALIZING
  START  WRITING ./tiledb-data/pbmc3k_processed
    START  WRITING ./tiledb-data/pbmc3k_processed/X/data
    FINISH WRITING ./tiledb-data/pbmc3k_processed/X/data
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
FINISH SCGroup.from_h5ad ./anndata/pbmc3k_processed.h5ad -> ./tiledb-data/pbmc3k_processed
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
