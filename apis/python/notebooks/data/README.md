How to regenerate:

* Source data -- download to local -- https://github.com/chanzuckerberg/cellxgene/raw/main/example-dataset/pbmc3k.h5ad
* Use `tiledbsoma.io.from_h5ad`
  * Remove previous `sparse/pbmc3k` and `dense/pbmc3k`
  * Comment out the call to the uns-writing logic (we don't need it)
  * Use default `X_kind` for `sparse/pbmc3k`
  * Use `X_kind=DenseNDArray` for `dense/pbmc3k`
