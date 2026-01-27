# Write R Objects to SOMA

Various helpers to write R objects to SOMA.

## Usage

``` r
# S3 method for class 'DataFrame'
write_soma(
  x,
  uri,
  soma_parent,
  df_index = NULL,
  index_column_names = "soma_joinid",
  ...,
  ingest_mode = "write",
  platform_config = NULL,
  tiledbsoma_ctx = NULL,
  context = NULL,
  relative = TRUE
)

# S3 method for class 'Hits'
write_soma(
  x,
  uri,
  soma_parent,
  sparse = TRUE,
  type = NULL,
  transpose = FALSE,
  ...,
  ingest_mode = "write",
  platform_config = NULL,
  tiledbsoma_ctx = NULL,
  context = NULL,
  relative = TRUE
)

# S3 method for class 'character'
write_soma(
  x,
  uri,
  soma_parent,
  ...,
  key = NULL,
  ingest_mode = "write",
  platform_config = NULL,
  tiledbsoma_ctx = NULL,
  context = NULL,
  relative = TRUE
)

# S3 method for class 'data.frame'
write_soma(
  x,
  uri,
  soma_parent,
  df_index = NULL,
  index_column_names = "soma_joinid",
  ...,
  key = NULL,
  ingest_mode = "write",
  platform_config = NULL,
  tiledbsoma_ctx = NULL,
  context = NULL,
  relative = TRUE
)

# S3 method for class 'IterableMatrix'
write_soma(
  x,
  uri,
  soma_parent,
  sparse = TRUE,
  type = NULL,
  transpose = FALSE,
  ...,
  key = NULL,
  ingest_mode = "write",
  shape = NULL,
  platform_config = NULL,
  tiledbsoma_ctx = NULL,
  context = NULL,
  relative = TRUE
)

# S3 method for class 'matrix'
write_soma(
  x,
  uri,
  soma_parent,
  sparse = TRUE,
  type = NULL,
  transpose = FALSE,
  ...,
  key = NULL,
  ingest_mode = "write",
  shape = NULL,
  platform_config = NULL,
  tiledbsoma_ctx = NULL,
  context = NULL,
  relative = TRUE
)

# S3 method for class 'Matrix'
write_soma(
  x,
  uri,
  soma_parent,
  sparse = TRUE,
  type = NULL,
  transpose = FALSE,
  ...,
  key = NULL,
  ingest_mode = "write",
  shape = NULL,
  platform_config = NULL,
  tiledbsoma_ctx = NULL,
  context = NULL,
  relative = TRUE
)

# S3 method for class 'TsparseMatrix'
write_soma(
  x,
  uri,
  soma_parent,
  type = NULL,
  transpose = FALSE,
  ...,
  key = NULL,
  ingest_mode = "write",
  shape = NULL,
  platform_config = NULL,
  tiledbsoma_ctx = NULL,
  context = NULL,
  relative = TRUE
)
```

## Arguments

- x:

  An object.

- uri:

  URI for resulting SOMA object.

- soma_parent:

  The parent [collection](SOMACollection.md) (eg. a
  [`SOMACollection`](SOMACollection.md),
  [`SOMAExperiment`](SOMAExperiment.md), or
  [`SOMAMeasurement`](SOMAMeasurement.md)).

- df_index:

  The name of the column in `x` with the index (row names); by default,
  will automatically add the row names of `x` to an attribute named
  “`index`” to the resulting [`SOMADataFrame`](SOMADataFrame.md).

- index_column_names:

  Names of columns in `x` to index in the resulting SOMA object.

- ...:

  Arguments passed to other methods

- ingest_mode:

  Ingestion mode when creating the SOMA; choose from:

  - “`write`”: create a new SOMA and error if it already exists.

  - “`resume`”: attempt to create a new SOMA; if it already exists,
    simply open it for writing.

- platform_config:

  Optional [platform configuration](PlatformConfig.md).

- tiledbsoma_ctx:

  Optional (DEPRECATED) [`SOMATileDBContext`](SOMATileDBContext.md).

- context:

  Optional `SOMAContext` object used for TileDB operations. If a context
  is not provided, then the default context will be used. Call
  `set_default_context` once before other SOMA operations to configure
  the default context.

- relative:

  **\[Internal use only\]** Is `uri` relative or absolute.

- sparse:

  Create a [sparse](SOMASparseNDArray.md) or
  [dense](SOMADenseNDArray.md) array from `x`.

- type:

  [Arrow type](https://arrow.apache.org/docs/r/reference/data-type.html)
  for encoding `x` (eg.
  [`arrow::int32()`](https://arrow.apache.org/docs/r/reference/data-type.html));
  by default, attempts to determine arrow type with
  [`arrow::infer_type()`](https://arrow.apache.org/docs/r/reference/infer_type.html).

- transpose:

  Transpose `x` before writing.

- key:

  Optionally register the resulting `SOMADataFrame` in `soma_parent` as
  `key`; pass `NULL` to prevent registration to handle manually.

- shape:

  A vector of two positive integers giving the on-disk shape of the
  array; defaults to `dim(x)`.

## Value

The resulting SOMA [array](SOMASparseNDArray.md) or [data
frame](SOMADataFrame.md), returned opened for write.

## Writing Character Arrays

[Characters](https://rdrr.io/r/base/character.html) are written out as
[`SOMADataFrame`](SOMADataFrame.md)`s` with one attribute: “`values`”;
additionally one bit of array-level metadata is added:

- “`soma_uns_outgest_hint`” with a value of “`array_1d`”.

## Writing Data Frames

[Data frames](https://rdrr.io/r/base/data.frame.html) are written out as
[`SOMADataFrame`](SOMADataFrame.md)`s`. The following transformations
are applied to `x`:

- row names are added to a column in `x` entitled “`index`”, “`_index`”,
  or a random name if either option is already present in `x`.

- a column “`soma_joinid`” will be automatically added going from
  `[0, nrow(x) - 1]` encoded as [64-bit
  integers](https://rdrr.io/pkg/bit64/man/bit64-package.html).

The array type for each column will be determined by
[`arrow::infer_type()`](https://arrow.apache.org/docs/r/reference/infer_type.html);
if any column contains a
[non-atomic](https://rdrr.io/r/base/is.recursive.html) type (excluding
[factors](https://rdrr.io/r/base/factor.html),
[`complex`](https://rdrr.io/r/base/complex.html)`es`,and
[`raw`](https://rdrr.io/r/base/raw.html)`s`), the code will error out.

## Writing Dense Matrices

Dense matrices are written as two-dimensional [dense
arrays](SOMADenseNDArray.md). The overall shape of the array is
determined by `dim(x)` and the type of the array is determined by `type`
or
[`arrow::infer_type`](https://arrow.apache.org/docs/r/reference/infer_type.html)`(x)`.

## Writing Sparse Matrices

Sparse matrices are written out as two-dimensional [TileDB sparse
arrays](SOMASparseNDArray.md) in [COO
format](https://en.wikipedia.org/wiki/Sparse_matrix#Coordinate_list_(COO)):

- the row indices (“`i`”) are written out as “`soma_dim_0`”.

- the column indices (“`j`”) are written out as “`soma_dim_1`”.

- the non-zero values (“`x`”) are written out as “`soma_data`”.

The array type is determined by `type`, or
[`arrow::infer_type`](https://arrow.apache.org/docs/r/reference/infer_type.html)`(slot(x, "x"))`.

## Examples

``` r
if (FALSE) { # requireNamespace("withr", quietly = TRUE) && requireNamespace("SeuratObject", quietly = TRUE) && requireNamespace("S4Vectors", quietly = TRUE)
# Write a Bioconductor S4 DataFrame object to a SOMA
uri <- withr::local_tempfile(pattern = "s4-data-frame")
data("pbmc_small", package = "SeuratObject")
obs <- suppressWarnings(SeuratObject::UpdateSeuratObject(pbmc_small))[[]]
head(obs <- as(obs, "DataFrame"))

(sdf <- write_soma(obs, uri, soma_parent = NULL, relative = FALSE))

sdf$close()
}
if (FALSE) { # requireNamespace("withr", quietly = TRUE) && requireNamespace("S4Vectors", quietly = TRUE)
# Write a Bioconductor SelfHits object to a SOMA
uri <- withr::local_tempfile(pattern = "hits")
(hits <- S4Vectors::SelfHits(
  c(2, 3, 3, 3, 3, 3, 4, 4, 4),
  c(4, 3, 2:4, 2, 2:3, 2),
  4,
  x = stats::rnorm(9L)
))

(arr <- write_soma(hits, uri, soma_parent = NULL, relative = FALSE))

arr$close()
}
# Write a character vector to a SOMA
uri <- withr::local_tempfile(pattern = "character")
(sdf <- write_soma(letters, uri, soma_parent = NULL, relative = FALSE))
#> <SOMADataFrame>
#>   uri: /tmp/Rtmpfr8mYm/character2a7e71bc2339
#>   dimensions: soma_joinid 
#>   attributes: values 

sdf$close()
# Write a data.frame to a SOMA
uri <- withr::local_tempfile(pattern = "data-frame")
data("pbmc_small", package = "SeuratObject")
head(obs <- suppressWarnings(SeuratObject::UpdateSeuratObject(pbmc_small))[[]])
#> Validating object structure
#> Updating object slots
#> Ensuring keys are in the proper structure
#> Updating matrix keys for DimReduc ‘pca’
#> Updating matrix keys for DimReduc ‘tsne’
#> Ensuring keys are in the proper structure
#> Ensuring feature names don't have underscores or pipes
#> Updating slots in RNA
#> Updating slots in RNA_snn
#> Setting default assay of RNA_snn to RNA
#> Updating slots in pca
#> Updating slots in tsne
#> Setting tsne DimReduc to global
#> Setting assay used for NormalizeData.RNA to RNA
#> Setting assay used for ScaleData.RNA to RNA
#> Setting assay used for RunPCA.RNA to RNA
#> Setting assay used for BuildSNN.RNA.pca to RNA
#> No assay information could be found for FindClusters
#> Setting assay used for RunTSNE.pca to RNA
#> Setting assay used for JackStraw.RNA.pca to RNA
#> Setting assay used for ScoreJackStraw.pca to RNA
#> Setting assay used for ProjectDim.RNA.pca to RNA
#> Setting assay used for FindVariableFeatures.RNA to RNA
#> Validating object structure for Assay ‘RNA’
#> Validating object structure for Graph ‘RNA_snn’
#> Validating object structure for DimReduc ‘pca’
#> Validating object structure for DimReduc ‘tsne’
#> Object representation is consistent with the most current Seurat version
#>                   orig.ident nCount_RNA nFeature_RNA RNA_snn_res.0.8
#> ATGCCAGAACGACT SeuratProject         70           47               0
#> CATGGCCTGTGCAT SeuratProject         85           52               0
#> GAACCTGATGAACC SeuratProject         87           50               1
#> TGACTGGATTCTCA SeuratProject        127           56               0
#> AGTCAGACTGCACA SeuratProject        173           53               0
#> TCTGATACACGTGT SeuratProject         70           48               0
#>                letter.idents groups RNA_snn_res.1
#> ATGCCAGAACGACT             A     g2             0
#> CATGGCCTGTGCAT             A     g1             0
#> GAACCTGATGAACC             B     g2             0
#> TGACTGGATTCTCA             A     g2             0
#> AGTCAGACTGCACA             A     g2             0
#> TCTGATACACGTGT             A     g1             0

(sdf <- write_soma(obs, uri, soma_parent = NULL, relative = FALSE))
#> <SOMADataFrame>
#>   uri: /tmp/Rtmpfr8mYm/data-frame2a7e43577f1d
#>   dimensions: soma_joinid 
#>   attributes: orig.ident, nCount_RNA, nFeature_RNA, RNA_snn_res.0.8, letter.idents, groups,... 

sdf$close()
if (FALSE) { # requireNamespace("withr", quietly = TRUE) && requireNamespace("SeuratObject", quietly = TRUE) && requireNamespace("BPCells", quietly = TRUE)
# Write a BPCells `IterableMatrix` to a SOMA
}
# Write a matrix to a SOMA
uri <- withr::local_tempfile(pattern = "matrix")
mat <- matrix(stats::rnorm(25L), nrow = 5L, ncol = 5L)
(arr <- write_soma(mat, uri, soma_parent = NULL, sparse = FALSE, relative = FALSE))
#> <SOMADenseNDArray>
#>   uri: /tmp/Rtmpfr8mYm/matrix2a7e24d892d0
#>   dimensions: soma_dim_0, soma_dim_1 
#>   attributes: soma_data 

arr$close()
# Write a dense S4 Matrix to a SOMA
uri <- withr::local_tempfile(pattern = "s4-matrix")
mat <- Matrix::Matrix(stats::rnorm(25L), nrow = 5L, ncol = 5L)
(arr <- write_soma(mat, uri, soma_parent = NULL, sparse = FALSE, relative = FALSE))
#> <SOMADenseNDArray>
#>   uri: /tmp/Rtmpfr8mYm/s4-matrix2a7ec3adaea
#>   dimensions: soma_dim_0, soma_dim_1 
#>   attributes: soma_data 

arr$close()
# Write a TsparseMatrix to a SOMA
uri <- withr::local_tempfile(pattern = "tsparse-matrix")
mat <- Matrix::rsparsematrix(5L, 5L, 0.3, repr = "T")
(arr <- write_soma(mat, uri, soma_parent = NULL, relative = FALSE))
#> <SOMASparseNDArray>
#>   uri: /tmp/Rtmpfr8mYm/tsparse-matrix2a7e9969814
#>   dimensions: soma_dim_0, soma_dim_1 
#>   attributes: soma_data 

arr$close()

# Write a CsparseMatrix to a SOMA
uri <- withr::local_tempfile(pattern = "csparse-matrix")
mat <- Matrix::rsparsematrix(5L, 5L, 0.3, repr = "C")
(arr <- write_soma(mat, uri, soma_parent = NULL, relative = FALSE))
#> <SOMASparseNDArray>
#>   uri: /tmp/Rtmpfr8mYm/csparse-matrix2a7e1e00a84
#>   dimensions: soma_dim_0, soma_dim_1 
#>   attributes: soma_data 

arr$close()

# Write an RsparseMatrix to a SOMA
uri <- withr::local_tempfile(pattern = "rsparse-matrix")
mat <- Matrix::rsparsematrix(5L, 5L, 0.3, repr = "R")
(arr <- write_soma(mat, uri, soma_parent = NULL, relative = FALSE))
#> <SOMASparseNDArray>
#>   uri: /tmp/Rtmpfr8mYm/rsparse-matrix2a7e691c15e9
#>   dimensions: soma_dim_0, soma_dim_1 
#>   attributes: soma_data 

arr$close()
```
