# pbmc-small

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

Note:

The R `anndata` package stores the distances adjacency matrix in a non-standard location, e.g:

```
x = ad.read_h5ad('/Users/aaronwolen/Downloads/pbmc-small.h5ad')
FutureWarning: Moving element from .uns['neighbors']['distances'] to .obsp['distances'].
```

But the Python `anndata` package (which we use for ingest) recognizes and fixes the problem.  The
`pbmc-small.h5ad` file included here has that corrected by having read and written back the file
using Python `anndata`.

# Various X storage

`pbmc3k_processed.h5ad` has dense X (`numpy.ndarray`). To produce CSR and CSC variants, to test dense/CSR/CSC
chunked writes, the following were done:

```
import anndata
import scipy

ann = anndata.read_h5ad('pbmc-small.h5ad')
ann.raw = None
ann.uns = {}
ann.obsm = None
ann.varm = None
ann.obsp = None
ann.write_h5ad('pbmc-small-x-dense.h5ad')

ann = anndata.read_h5ad('pbmc-small.h5ad')
ann.raw = None
ann.uns = {}
ann.obsm = None
ann.varm = None
ann.obsp = None
ann.X = scipy.sparse.csr_matrix(ann.X)
ann.write_h5ad('pbmc-small-x-csr.h5ad')

ann = anndata.read_h5ad('pbmc-small.h5ad')
ann.raw = None
ann.uns = {}
ann.obsm = None
ann.varm = None
ann.obsp = None
ann.X = scipy.sparse.csc_matrix(ann.X)
ann.write_h5ad('pbmc-small-x-csc.h5ad')
```
