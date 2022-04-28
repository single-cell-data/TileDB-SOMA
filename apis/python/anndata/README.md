Provenance of `pbmc-small.md` is from the Seurat R package:

```
> remotes::install_github("mojaveazure/seurat-disk")
> library(SeuratDisk)
> library(anndata)
>
> data("pbmc_small", package = "SeuratObject")
> SaveH5Seurat(pbmc_small, filename='pbmc-small.h5Seurat')
> Convert("pbmc-small.h5Seurat", dest="h5ad")
```
