# Write a [`Seurat`](https://satijalab.github.io/seurat-object/reference/Seurat-class.html) object to a SOMA

Write a
[`Seurat`](https://satijalab.github.io/seurat-object/reference/Seurat-class.html)
object to a SOMA

## Usage

``` r
# S3 method for class 'Seurat'
write_soma(
  x,
  uri,
  ...,
  ingest_mode = "write",
  platform_config = NULL,
  tiledbsoma_ctx = NULL,
  context = NULL
)
```

## Arguments

- x:

  A
  [`Seurat`](https://satijalab.github.io/seurat-object/reference/Seurat-class.html)
  object.

- uri:

  URI for resulting SOMA object.

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

## Value

The URI to the resulting [`SOMAExperiment`](SOMAExperiment.md) generated
from the data contained in `x`.

## Writing Cell-Level Metadata

Cell-level metadata is written out as a [data frame](SOMADataFrame.md)
called “`obs`” at the [`experiment`](SOMAExperiment.md) level.

## Writing v3 [`Assay`](https://satijalab.github.io/seurat-object/reference/Assay-class.html)`s`

Seurat
[`Assay`](https://satijalab.github.io/seurat-object/reference/Assay-class.html)
objects are written out as individual
[measurements](SOMAMeasurement.md):

- the “`data`” matrix is written out as a [sparse
  array](SOMASparseNDArray.md) called “`data`” within the “`X`” group.

- the “`counts`” matrix, if not
  [empty](https://satijalab.github.io/seurat-object/reference/IsMatrixEmpty.html),
  is written out as a [sparse array](SOMASparseNDArray.md) called
  “`counts`” within the “`X`” group.

- the “`scale.data`” matrix, if not
  [empty](https://satijalab.github.io/seurat-object/reference/IsMatrixEmpty.html),
  is written out as a [sparse array](SOMASparseNDArray.md) called
  “`scale_data`” within the “`X`” group.

- feature-level metadata is written out as a [data
  frame](SOMADataFrame.md) called “`var`”.

Expression matrices are transposed (cells as rows) prior to writing. All
other slots, including results from extended assays (eg. `SCTAssay`,
`ChromatinAssay`) are lost.

### Performance Considerations

Ingestion of very large dense layers, such as `scale.data`, can be
memory intensive. For better performance, users can remove these layers
prior to ingestion and regenerate them after export, or ingest them
separately as dense arrays for those who need to persist the exact
matrix

    # Using SeuratObject v5 syntax on a v3 `Assay`
    # Cache the layer for separate ingestion, skip if planning to regenerate
    mat <- object[["ASSAY"]]$scale.data

    # Remove the `scale.data` layer
    object[["ASSAY"]]$scale.data <- NULL

    # Ingest the smaller object
    uri <- write_soma(object, "/path/to/soma")

    # Ingest the `scale.data` layer densely; needed only if persistence
    # of the data is paramount
    # Pad the `scale.data` layer so that its soma join IDs match the experiment
    padded <- matrix(
      data = vector("numeric", length = prod(dim(object[["ASSAY"]]))),
      nrow = nrow(object[["ASSAY"]]),
      ncol = ncol(object[["ASSAY"]])
    )
    rowidx <- match(rownames(mat), rownames(object[["ASSAY"]]))
    colidx <- match(colnames(mat), colnames(object[["ASSAY"]]))
    padded[rowidx, colidx] <- mat

    # Use `write_soma()` to ingest densely and register it within the `uns`
    # collection; this may need to be created manually if the original
    # object does not contain command logs
    exp <- SOMAExperimentOpen(uri, "WRITE")
    if (!match("uns", exp$names(), nomatch = 0L)) {
      # For `tiledb://` URIs, set the URI for the new collection manually rather
      # than relying on `file.path()`
      uns <- SOMACollectionCreate(file.path(exp$uri, "uns"))
      exp$add_new_collection(uns, "uns")
    }
    arr <- write_soma(
      padded,
      "scale_data",
      soma_parent = exp$get("uns"),
      sparse = FALSE,
      key = "scale_data"
    )
    arr$close()
    exp$close()

Please note that dense arrays cannot be read in using the
[`SOMAExperimentAxisQuery`](SOMAExperimentAxisQuery.md) mechanism; use
[`SOMADenseNDArray`](SOMADenseNDArray.md)`$read_dense_matrix`,
remembering to transpose before adding back to a `Seurat` object

## Writing v5 `Assays`

Seurat v5
[`Assays`](https://satijalab.github.io/seurat-object/reference/Assay5-class.html)`s`
are written out as individual [measurements](SOMAMeasurement.md):

- the layer matrices are written out as [sparse
  arrays](SOMASparseNDArray.md) within the “`X`” group.

- feature-level metadata is written out as a [data
  frame](SOMADataFrame.md) called “`var`”.

Expression matrices are transposed (cells as rows) prior to writing. All
other slots, including results from extended assays (eg. `SCTAssay`,
`ChromatinAssay`) are lost.  
The following bits of metadata are written in various parts of the
measurement

- “`soma_ecosystem_seurat_assay_version`”: written at the measurement
  level; indicates the Seurat assay version. Set to “`v5`”.

- “`soma_ecosystem_seurat_v5_default_layers`”: written at the “`X`”
  group level; indicates the [default
  layers](https://satijalab.github.io/seurat-object/reference/DefaultLayer.html).

- “`soma_ecosystem_seurat_v5_ragged`”: written at the “`X/<layer>`”
  array level; with a value of “`ragged`”, indicates whether or not the
  layer is ragged.

- “`soma_r_type_hint`”: written at the “`X/<layer>`” array level;
  indicates the R class and defining package (for S4 classes) of the
  original layer.

## Writing [`DimReduc`](https://satijalab.github.io/seurat-object/reference/DimReduc-class.html)`s`

Seurat
[`DimReduc`](https://satijalab.github.io/seurat-object/reference/DimReduc-class.html)
objects are written out to the “`obsm`” and “`varm`” groups of a
[measurement](SOMAMeasurement.md):

- cell embeddings are written out as a [sparse
  matrix](SOMASparseNDArray.md) in the “`obsm`” group.

- feature loadings, if not
  [empty](https://satijalab.github.io/seurat-object/reference/IsMatrixEmpty.html),
  are written out as a [sparse matrix](SOMASparseNDArray.md) in the
  “`varm`” groups; loadings are padded with `NAs` to include all
  features.

Dimensional reduction names are translated to AnnData-style names (eg.
“`pca`” becomes `X_pca` for embeddings and “`PCs`” for loadings). All
other slots, including projected feature loadings and jackstraw
information, are lost.

## Writing [`Graph`](https://satijalab.github.io/seurat-object/reference/Graph-class.html)`s`

Seurat
[`Graph`](https://satijalab.github.io/seurat-object/reference/Graph-class.html)
objects are written out as [sparse matrices](SOMASparseNDArray.md) to
the “`obsp`” group of a [measurement](SOMAMeasurement.md).

## Writing [`SeuratCommand`](https://satijalab.github.io/seurat-object/reference/SeuratCommand-class.html)`s`

Seurat [command
logs](https://satijalab.github.io/seurat-object/reference/SeuratCommand-class.html)
are written out as [data frames](SOMADataFrame.md) to the
“`seurat_commands`” group of a [collection](SOMACollection.md).

## Examples

``` r
# \donttest{
uri <- withr::local_tempfile(pattern = "pbmc-small")

data("pbmc_small", package = "SeuratObject")
suppressWarnings(pbmc_small <- SeuratObject::UpdateSeuratObject(pbmc_small))
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

uri <- write_soma(pbmc_small, uri)

(exp <- SOMAExperimentOpen(uri))
#> <SOMAExperiment>
#>   uri: /tmp/Rtmpfr8mYm/pbmc-small2a7e32b82591
exp$obs
#> <SOMADataFrame>
#>   uri: file:///tmp/Rtmpfr8mYm/pbmc-small2a7e32b82591/obs
#>   dimensions: soma_joinid 
#>   attributes: orig.ident, nCount_RNA, nFeature_RNA, RNA_snn_res.0.8, letter.idents, groups,... 
exp$get("uns")$get("seurat_commands")$names()
#>  [1] "BuildSNN.RNA.pca"         "FindClusters"            
#>  [3] "FindVariableFeatures.RNA" "JackStraw.RNA.pca"       
#>  [5] "NormalizeData.RNA"        "ProjectDim.RNA.pca"      
#>  [7] "RunPCA.RNA"               "RunTSNE.pca"             
#>  [9] "ScaleData.RNA"            "ScoreJackStraw.pca"      
(ms <- exp$ms$get("RNA"))
#> <SOMAMeasurement>
#>   uri: file:///tmp/Rtmpfr8mYm/pbmc-small2a7e32b82591/ms/RNA
ms$var
#> <SOMADataFrame>
#>   uri: file:///tmp/Rtmpfr8mYm/pbmc-small2a7e32b82591/ms/RNA/var
#>   dimensions: soma_joinid 
#>   attributes: vst.mean, vst.variance, vst.variance.expected, vst.variance.standardized, vst... 
ms$X$names()
#> [1] "counts"     "data"       "scale_data"
ms$obsm$names()
#> [1] "X_pca"  "X_tsne"
ms$varm$names()
#> [1] "PCs"
ms$obsp$names()
#> [1] "RNA_snn"

exp$close()
# }
```
