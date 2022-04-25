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

# Discussion

* This is currently on the `python-write-path-v2` branch of my fork repo [https://github.com/johnkerl/TileDB-SingleCell/tree/python-write-path-v2](https://github.com/johnkerl/TileDB-SingleCell/tree/python-write-path-v2)
* `anndata/pbmc3k_processed.h5ad` is a file CZI recommended to us
* The v1 code `v1/anndata_to_tiledb.py` does not successfully process the CZI-recommended `pbmc3k_processed.h5ad`
* In the HD5 file we read, there is `X`, `obs`, `obsm`, `raw.X`, `raw.var`, `uns`, `var`, `varm`
* In the TileDB group we write, we want to have all of those.
