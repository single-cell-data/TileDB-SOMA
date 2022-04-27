# Overview

Compare `pbmc_small` as written by

* `tiledbsc-r`: [https://github.com/TileDB-Inc/tiledbsc](https://github.com/TileDB-Inc/tiledbsc)
* `tiledbsc-py`: [https://github.com/single-cell-data/TileDB-SingleCell/tree/main/apis/python](https://github.com/single-cell-data/TileDB-SingleCell/tree/main/apis/python)

# Setup

## R setup

```
install.packages("Seurat")
remotes::install_github("mojaveazure/seurat-disk")
install.packages('anndata')
```

## Python setup

```
pip install anndata
```

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

## Notes

* _Unsure if the R `anndata` package is a route we should be using_

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
* This is because the rows/columns are transposed in R (seurat/bioconductor) vs. Python (anndata). This is not unexpected but how/whether we should address it is unclear.
* In R the data matrices are genes by samples (or cells). Whereas in Python/anndata they are samples by genes.

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

## Matrix-data level

### X/data

```
> arr <- tiledb_array('pbmc-small-tiledbsc-r/scgroup_RNA/X/data', return_as='data.frame')
> df <- arr[]
> head(df$value)
[1] 4.615121 4.208048 3.976987 3.997759 3.485709 4.753095
> head(df$obs_id)
[1] "AATGTTGACAGTCA" "ACAGGTACTGGTGT" "AGAGATGATCTCGC" "AGATATACCCGTAA"
[5] "CATTACACCAACTG" "GAACCTGATGAACC"
> head(df$var_id)
[1] "ACAP1" "ACAP1" "ACAP1" "ACAP1" "ACAP1" "ACAP1"
```

```
>>> import tiledb
>>> arr = tiledb.open('pbmc-small-tiledbsc-r/scgroup_RNA/X/data')
>>> df = arr[]
>>> df.keys()
odict_keys(['value', 'var_id', 'obs_id'])
>>> df['value']
array([4.61512052, 4.20804766, 3.97698683, ..., 4.2412282 , 4.17596394, 4.37877332])
>>> df['obs_id']
array([b'AATGTTGACAGTCA', b'ACAGGTACTGGTGT', b'AGAGATGATCTCGC', ..., b'CTGCCAACAGGAGC', b'TAGGGACTGAACTC', b'TGACTGGATTCTCA'], dtype=object)
>>> df['var_id']
array([b'ACAP1', b'ACAP1', b'ACAP1', ..., b'ZNF76', b'ZNF76', b'ZNF76'], dtype=object)
```

```
> arr <- tiledb_array('pbmc-small-tiledbsc-py/X/data', return_as='data.frame')
> df <- arr[]
> head(df$data)
[1] -0.3730316  2.5200651 -0.3730316 -0.3730316 -0.3730316  2.6456511
> head(df$obs_id)
[1] "AAATTCGAATCACG" "AAATTCGAATCACG" "AAATTCGAATCACG" "AAATTCGAATCACG"
[5] "AAATTCGAATCACG" "AAATTCGAATCACG"
> head(df$var_id)
[1] "AKR1C3"   "CA2"      "CD1C"     "GNLY"     "HLA-DPB1" "HLA-DQA1"
```

```
>>> import tiledb
>>> arr = tiledb.open('pbmc-small-tiledbsc-py/X/data')
>>> df = arr[:]
>>> df.keys()
odict_keys(['data', 'obs_id', 'var_id'])
>>> df['data']
array([-0.37303159,  2.52006507, -0.37303159, ...,  0.45850021, -1.039994  , -1.039994  ])
>>> df['obs_id']
array([b'AAATTCGAATCACG', b'AAATTCGAATCACG', b'AAATTCGAATCACG', ..., b'TTTAGCTGTACTCT', b'TTTAGCTGTACTCT', b'TTTAGCTGTACTCT'], dtype=object)
>>> df['var_id']
array([b'AKR1C3', b'CA2', b'CD1C', ..., b'TREML1', b'TUBB1', b'VDAC3'], dtype=object)
```

### obs

```
> arr <- tiledb_array('pbmc-small-tiledbsc-r/scgroup_RNA/obs', return_as='data.frame')
> df = arr[]
> df$
df$RNA_snn_res.0.8  df$groups           df$nCount_RNA       df$obs_id
df$RNA_snn_res.1    df$letter.idents    df$nFeature_RNA     df$orig.ident
>
> head(df$obs_id)
[1] "AAATTCGAATCACG" "AAGCAAGAGCTTAG" "AAGCGACTTTGACG" "AATGCGTGGACGGA"
[5] "AATGTTGACAGTCA" "ACAGGTACTGGTGT"
> head(df$`RNA_snn_res.0.8`)
[1] "1" "0" "1" "1" "0" "0"
> head(df$letter.idents)
[1] "B" "A" "B" "B" "A" "A"
```

```
> arr <- tiledb_array('pbmc-small-tiledbsc-py/obs', return_as='data.frame')
> df = arr[]
> df$
df$RNA_snn_res.0.8  df$__tiledb_rows    df$letter.idents    df$nFeature_RNA
df$RNA_snn_res.1    df$groups           df$nCount_RNA       df$orig.ident
>
> head(df$`__tiledb_rows`)
[1] "AAATTCGAATCACG" "AAGCAAGAGCTTAG" "AAGCGACTTTGACG" "AATGCGTGGACGGA"
[5] "AATGTTGACAGTCA" "ACAGGTACTGGTGT"
> head(df$`RNA_snn_res.0.8`)
[1] 1 0 1 1 0 0
> head(df$`letter.idents`)
[1] 1 0 1 1 0 0
```
