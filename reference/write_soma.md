# Write a SOMA Object from an R Object

Convert R objects to their appropriate SOMA counterpart function and
methods can be written for it to provide a high-level R \\\rightarrow\\
SOMA interface.

## Usage

``` r
write_soma(
  x,
  uri,
  ...,
  platform_config = NULL,
  tiledbsoma_ctx = NULL,
  context = NULL
)
```

## Arguments

- x:

  An object.

- uri:

  URI for resulting SOMA object.

- ...:

  Arguments passed to other methods

- platform_config:

  Optional [platform configuration](PlatformConfig.md).

- tiledbsoma_ctx:

  Optional (DEPRECATED) [`SOMATileDBContext`](SOMATileDBContext.md).

- context:

  Optional `SOMAContext` object used for TileDB operations. If a context
  is not provided, then the default context will be used. Call
  `set_default_context` once before other SOMA operations to configure
  the default context.

## Value

The URI to the resulting [`SOMAExperiment`](SOMAExperiment.md) generated
from the data contained in `x`.

## Known methods

- [Writing Seurat objects](write_soma.Seurat.md).

- [Writing SummarizedExperiment
  objects](write_soma.SummarizedExperiment.md).

- [Writing SingleCellExperiment
  objects](write_soma.SingleCellExperiment.md).

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
#>   uri: /tmp/Rtmpfr8mYm/character2a7e1a7eb7b1
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
#>   uri: /tmp/Rtmpfr8mYm/data-frame2a7e463f18f7
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
#>   uri: /tmp/Rtmpfr8mYm/matrix2a7e5d0777b4
#>   dimensions: soma_dim_0, soma_dim_1 
#>   attributes: soma_data 

arr$close()
# Write a dense S4 Matrix to a SOMA
uri <- withr::local_tempfile(pattern = "s4-matrix")
mat <- Matrix::Matrix(stats::rnorm(25L), nrow = 5L, ncol = 5L)
(arr <- write_soma(mat, uri, soma_parent = NULL, sparse = FALSE, relative = FALSE))
#> <SOMADenseNDArray>
#>   uri: /tmp/Rtmpfr8mYm/s4-matrix2a7e5ce13aff
#>   dimensions: soma_dim_0, soma_dim_1 
#>   attributes: soma_data 

arr$close()
# Write a TsparseMatrix to a SOMA
uri <- withr::local_tempfile(pattern = "tsparse-matrix")
mat <- Matrix::rsparsematrix(5L, 5L, 0.3, repr = "T")
(arr <- write_soma(mat, uri, soma_parent = NULL, relative = FALSE))
#> <SOMASparseNDArray>
#>   uri: /tmp/Rtmpfr8mYm/tsparse-matrix2a7e3844124e
#>   dimensions: soma_dim_0, soma_dim_1 
#>   attributes: soma_data 

arr$close()

# Write a CsparseMatrix to a SOMA
uri <- withr::local_tempfile(pattern = "csparse-matrix")
mat <- Matrix::rsparsematrix(5L, 5L, 0.3, repr = "C")
(arr <- write_soma(mat, uri, soma_parent = NULL, relative = FALSE))
#> <SOMASparseNDArray>
#>   uri: /tmp/Rtmpfr8mYm/csparse-matrix2a7e7b919c5b
#>   dimensions: soma_dim_0, soma_dim_1 
#>   attributes: soma_data 

arr$close()

# Write an RsparseMatrix to a SOMA
uri <- withr::local_tempfile(pattern = "rsparse-matrix")
mat <- Matrix::rsparsematrix(5L, 5L, 0.3, repr = "R")
(arr <- write_soma(mat, uri, soma_parent = NULL, relative = FALSE))
#> <SOMASparseNDArray>
#>   uri: /tmp/Rtmpfr8mYm/rsparse-matrix2a7e44ce477f
#>   dimensions: soma_dim_0, soma_dim_1 
#>   attributes: soma_data 

arr$close()
```
