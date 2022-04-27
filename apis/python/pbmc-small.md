# Setup

## R setup

```
install.packages("Seurat")
remotes::install_github("mojaveazure/seurat-disk")
install.packages('anndata')
```

`tiledbsc-r`: [https://github.com/TileDB-Inc/tiledbsc](https://github.com/TileDB-Inc/tiledbsc)

## Python setup

```
pip install anndata
```

`tiledbsc-py`: [https://github.com/single-cell-data/TileDB-SingleCell/tree/main/apis/python](https://github.com/single-cell-data/TileDB-SingleCell/tree/main/apis/python)

# Write pbmc_small datasets

## Convert Seurat to tiledbsc-r

```
library(tiledbsc)

data("pbmc_small", package = "SeuratObject")
scdataset <- SCDataset$new(uri = "pbmc-small-tiledbsc-r")
scdataset$from_seurat(object = pbmc_small)
```

## Convert Seurat to anndata

```
library(anndata)
library(SeuratDisk)

data("pbmc_small", package = "SeuratObject")
SaveH5Seurat(pbmc_small, filename='pbmc-small.h5Seurat')
Convert("pbmc-small.h5Seurat", dest="pbmc-small.h5ad")
```

## Convert anndata to tiledbsc-py

Or use `ingestory.py` wrapper script which does this.

```
import tiledbsc

scdataset = tiledbsc.SCGroup("pbmc-small-tiledbsc-py", verbose=True)
scdataset.from_h5ad("pbmc-small.h5ad")
```

# Compare pbmc_small datasets

## Group-tree level

```
import tiledb

print(tiledb.group.Group('pbmc-small-tiledbsc-r/scgroup_RNA')._dump(True))

print(tiledb.group.Group('pbmc-small-tiledbsc-py')._dump(True))
```

```
>>> import tiledb
>>>
>>> print(tiledb.group.Group('pbmc-small-tiledbsc-r/scgroup_RNA')._dump(True))
scgroup_RNA GROUP
|-- X GROUP
|------ data ARRAY
|------ counts ARRAY
|------ scale.data ARRAY
|-- obs ARRAY
|-- var ARRAY
|-- obsm GROUP
|------ dimreduction_pca ARRAY
|------ dimreduction_tsne ARRAY
|-- varm GROUP
|------ dimreduction_pca ARRAY
|-- varp GROUP
|-- obsp GROUP
|------ graph_snn ARRAY
|-- misc GROUP

>>>
>>>
>>> print(tiledb.group.Group('pbmc-small-tiledbsc-py')._dump(True))
pbmc-small-tiledbsc-py GROUP
|-- X GROUP
|------ data ARRAY
|------ raw ARRAY
|-- obs ARRAY
|-- var ARRAY
|-- obsm GROUP
|------ X_pca ARRAY
|------ X_tsne ARRAY
|-- varm GROUP
|------ PCs ARRAY
|-- obsp GROUP
|------ distances ARRAY
```

## Array-schema level

First note: display TileDB array schemas using `print(tiledb.open(uri).schema)`. Command-line
keystroke-saver `show-array-schema.py`:

```
$ cat $(show-array-schema.py)
#!/usr/bin/env python

import sys, os
import tiledb

if len(sys.argv) < 2:
    print("Usage: %s {one or more array URIs}" % (sys.argv[0]), file=sys.stderr)
    sys.exit(1)

for uri in sys.argv[1:]:
    print(tiledb.open(uri).schema)
```

### X/data

```
$ show-array-schema.py pbmc-small-tiledbsc-r/scgroup_RNA/X/data/
ArraySchema(
  domain=Domain(*[
    Dim(name='var_id', domain=(None, None), tile=None, dtype='|S0', var=True),
    Dim(name='obs_id', domain=(None, None), tile=None, dtype='|S0', var=True),
  ]),
  attrs=[
    Attr(name='value', dtype='float64', var=False, nullable=False, filters=FilterList([ZstdFilter(level=-1), ])),
  ],
  cell_order='row-major',
  tile_order='row-major',
  capacity=10000,
  sparse=True,
  allows_duplicates=True,
)
```

```
$ show-array-schema.py pbmc-small-tiledbsc-py/X/data
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
  allows_duplicates=True,
)
```

**Note:**

* `obs_id` and `var_id` are reversed between R and Python

### obs

```
$ show-array-schema.py pbmc-small-tiledbsc-r/scgroup_RNA/obs
ArraySchema(
  domain=Domain(*[
    Dim(name='obs_id', domain=(None, None), tile=None, dtype='|S0', var=True),
  ]),
  attrs=[
    Attr(name='orig.ident', dtype='ascii', var=True, nullable=False, filters=FilterList([ZstdFilter(level=-1), ])),
    Attr(name='nCount_RNA', dtype='float64', var=False, nullable=False, filters=FilterList([ZstdFilter(level=-1), ])),
    Attr(name='nFeature_RNA', dtype='int32', var=False, nullable=False, filters=FilterList([ZstdFilter(level=-1), ])),
    Attr(name='RNA_snn_res.0.8', dtype='ascii', var=True, nullable=False, filters=FilterList([ZstdFilter(level=-1), ])),
    Attr(name='letter.idents', dtype='ascii', var=True, nullable=False, filters=FilterList([ZstdFilter(level=-1), ])),
    Attr(name='groups', dtype='ascii', var=True, nullable=False, filters=FilterList([ZstdFilter(level=-1), ])),
    Attr(name='RNA_snn_res.1', dtype='ascii', var=True, nullable=False, filters=FilterList([ZstdFilter(level=-1), ])),
  ],
  cell_order='row-major',
  tile_order='row-major',
  capacity=10000,
  sparse=True,
  allows_duplicates=True,
)
```

```
$ show-array-schema.py pbmc-small-tiledbsc-py/obs
ArraySchema(
  domain=Domain(*[
    Dim(name='__tiledb_rows', domain=(None, None), tile=None, dtype='|S0', var=True, filters=FilterList([ZstdFilter(level=-1), ])),
  ]),
  attrs=[
    Attr(name='orig.ident', dtype='int32', var=False, nullable=False, filters=FilterList([ZstdFilter(level=-1), ])),
    Attr(name='nCount_RNA', dtype='float64', var=False, nullable=False, filters=FilterList([ZstdFilter(level=-1), ])),
    Attr(name='nFeature_RNA', dtype='int32', var=False, nullable=False, filters=FilterList([ZstdFilter(level=-1), ])),
    Attr(name='RNA_snn_res.0.8', dtype='int32', var=False, nullable=False, filters=FilterList([ZstdFilter(level=-1), ])),
    Attr(name='letter.idents', dtype='int32', var=False, nullable=False, filters=FilterList([ZstdFilter(level=-1), ])),
    Attr(name='groups', dtype='<U0', var=True, nullable=False, filters=FilterList([ZstdFilter(level=-1), ])),
    Attr(name='RNA_snn_res.1', dtype='int32', var=False, nullable=False, filters=FilterList([ZstdFilter(level=-1), ])),
  ],
  cell_order='row-major',
  tile_order='row-major',
  capacity=10000,
  sparse=True,
  allows_duplicates=False,
)
```

### var

```
$ show-array-schema.py pbmc-small-tiledbsc-r/scgroup_RNA/var
ArraySchema(
  domain=Domain(*[
    Dim(name='var_id', domain=(None, None), tile=None, dtype='|S0', var=True),
  ]),
  attrs=[
    Attr(name='vst.mean', dtype='float64', var=False, nullable=False, filters=FilterList([ZstdFilter(level=-1), ])),
    Attr(name='vst.variance', dtype='float64', var=False, nullable=False, filters=FilterList([ZstdFilter(level=-1), ])),
    Attr(name='vst.variance.expected', dtype='float64', var=False, nullable=False, filters=FilterList([ZstdFilter(level=-1), ])),
    Attr(name='vst.variance.standardized', dtype='float64', var=False, nullable=False, filters=FilterList([ZstdFilter(level=-1), ])),
    Attr(name='vst.variable', dtype='int32', var=False, nullable=False, filters=FilterList([ZstdFilter(level=-1), ])),
    Attr(name='highly_variable', dtype='int32', var=False, nullable=False, filters=FilterList([ZstdFilter(level=-1), ])),
  ],
  cell_order='row-major',
  tile_order='row-major',
  capacity=10000,
  sparse=True,
  allows_duplicates=True,
)
```

```
$ show-array-schema.py pbmc-small-tiledbsc-py/var
ArraySchema(
  domain=Domain(*[
    Dim(name='__tiledb_rows', domain=(None, None), tile=None, dtype='|S0', var=True, filters=FilterList([ZstdFilter(level=-1), ])),
  ]),
  attrs=[
    Attr(name='vst.mean', dtype='float64', var=False, nullable=False, filters=FilterList([ZstdFilter(level=-1), ])),
    Attr(name='vst.variance', dtype='float64', var=False, nullable=False, filters=FilterList([ZstdFilter(level=-1), ])),
    Attr(name='vst.variance.expected', dtype='float64', var=False, nullable=False, filters=FilterList([ZstdFilter(level=-1), ])),
    Attr(name='vst.variance.standardized', dtype='float64', var=False, nullable=False, filters=FilterList([ZstdFilter(level=-1), ])),
    Attr(name='vst.variable', dtype='int32', var=False, nullable=False, filters=FilterList([ZstdFilter(level=-1), ])),
  ],
  cell_order='row-major',
  tile_order='row-major',
  capacity=10000,
  sparse=True,
  allows_duplicates=False,
)
```

### obsm

```
$ show-array-schema.py pbmc-small-tiledbsc-r/scgroup_RNA/obsm/dimreduction_tsne
ArraySchema(
  domain=Domain(*[
    Dim(name='obs_id', domain=(None, None), tile=None, dtype='|S0', var=True),
  ]),
  attrs=[
    Attr(name='tSNE_1', dtype='float64', var=False, nullable=False, filters=FilterList([ZstdFilter(level=-1), ])),
    Attr(name='tSNE_2', dtype='float64', var=False, nullable=False, filters=FilterList([ZstdFilter(level=-1), ])),
  ],
  cell_order='row-major',
  tile_order='row-major',
  capacity=10000,
  sparse=True,
  allows_duplicates=True,
)
```

```
$ show-array-schema.py pbmc-small-tiledbsc-py/obsm/X_tsne/
ArraySchema(
  domain=Domain(*[
    Dim(name='__dim_0', domain=(0, 79), tile=80, dtype='uint64'),
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
```

# Notes

* Unsure if the R `anndata` package is a route we should be using
