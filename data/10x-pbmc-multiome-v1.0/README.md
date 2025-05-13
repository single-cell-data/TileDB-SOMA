# Data info

Small sample dataset for testing.

- source: 10x Genomics
  - https://www.10xgenomics.com/resources/datasets/10-k-human-pbm-cs-multiome-v-1-0-chromium-controller-1-standard-2-0-0
- license: Creative Commons By-Attribution
- subsetting from original, with scanpy:

```
import scanpy

ds = scanpy.read_10x_h5(
    "10k_PBMC_Multiome_nextgem_Chromium_Controller_raw_feature_bc_matrix.h5"
)
ds[:100, :100].write_h5ad("subset_1.h5ad")
```
