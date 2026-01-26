# Write a `SummarizedExperiment` object to a SOMA

Write a `SummarizedExperiment` object to a SOMA

## Usage

``` r
# S3 method for class 'SummarizedExperiment'
write_soma(
  x,
  uri,
  ms_name,
  ...,
  ingest_mode = "write",
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

- ms_name:

  Name for resulting measurement.

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

## Writing `colData`

`colData` is written out as a [data frame](SOMADataFrame.md) called
“`obs`” at the [`experiment`](SOMAExperiment.md) level.

## Writing Assay Matrices

Each assay matrix is written out as a [sparse
matrix](SOMASparseNDArray.md) within the `X` group of
[`measurement`](SOMAMeasurement.md) named `ms_name`. Names for assay
matrices within `X` are taken from the assay names. Assay matrices are
transposed (samples as rows) prior to writing.

## Writing `rowData`

`rowData` is written out as a [data frame](SOMADataFrame.md) called
“`var`” at the [`measurement`](SOMAMeasurement.md) level.

## Examples

``` r
if (FALSE) { # requireNamespace("withr", quietly = TRUE) && requireNamespace("SummarizedExperiment", quietly = TRUE)
# \donttest{
uri <- withr::local_tempfile(pattern = "summarized-experiment")

mat <- abs(Matrix::rsparsematrix(
  230L,
  80L,
  0.3,
  dimnames = list(paste0("feature_", seq_len(230)), paste0("cell_", seq_len(80)))
))
(se <- SummarizedExperiment::SummarizedExperiment(list(counts = mat, logcounts = log2(mat + 1L))))

uri <- write_soma(se, uri, ms_name = "RNA")

(exp <- SOMAExperimentOpen(uri))
exp$obs
(ms <- exp$ms$get("RNA"))
ms$var
ms$X$names()

exp$close()
# }
}
```
