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

## Notes

* Unsure if the R `anndata` package is a route we should be using
