## Creating data examples

To create the SOMAExperiments in `./data` from a python session do the following.

### Download 10X PBMC data from CELLxGENE test files

```bash
wget https://github.com/chanzuckerberg/cellxgene/raw/main/example-dataset/pbmc3k.h5ad
rm -r ./data/sparse/pbmc3k ./data/dense/pbmc3k
```

### Make `SOMAExperiment`

From python

```python
import anndata
import tiledbsoma.io
from tiledbsoma import SparseNDArray, DenseNDArray

adata = anndata.read("pbmc3k.h5ad")
del adata.uns

tiledbsoma.io.from_anndata("./data/dense/pbmc3k", adata, measurement_name = "RNA", X_kind = DenseNDArray)
tiledbsoma.io.from_anndata("./data/sparse/pbmc3k", adata, measurement_name = "RNA", X_kind = SparseNDArray)
"

```
