# Open a SOMA Object

Utility function to open the corresponding SOMA object given a URI
(lifecycle: maturing).

## Usage

``` r
SOMAOpen(
  uri,
  mode = "READ",
  platform_config = NULL,
  tiledbsoma_ctx = NULL,
  tiledb_timestamp = NULL
)
```

## Arguments

- uri:

  URI for the TileDB object.

- mode:

  One of “`READ`” or “`WRITE`”

- platform_config:

  Optional platform configuration.

- tiledbsoma_ctx:

  Optional SOMATileDBContext.

- tiledb_timestamp:

  Optional Datetime (POSIXct) for TileDB timestamp; defaults to the
  current time. If not `NULL`, all members accessed through the
  collection inherit the timestamp.

## Value

A SOMA object

## Examples

``` r
dir <- withr::local_tempfile(pattern = "soma-open")
dir.create(dir, recursive = TRUE)

uri <- extract_dataset("soma-exp-pbmc-small", dir)
(exp <- SOMAOpen(uri))
#> <SOMAExperiment>
#>   uri: /tmp/Rtmpf7VWBj/soma-open2b795e140fb/soma-exp-pbmc-small


uri <- extract_dataset("soma-dataframe-pbmc3k-processed-obs", dir)
(obs <- SOMAOpen(uri))
#> <SOMADataFrame>
#>   uri: /tmp/Rtmpf7VWBj/soma-open2b795e140fb/soma-dataframe-pbmc3k-processed-obs
#>   dimensions: soma_joinid 
#>   attributes: orig.ident, nCount_RNA, nFeature_RNA, seurat_annotations, percent.mt, RNA_snn... 
```
