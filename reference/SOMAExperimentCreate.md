# Create a SOMA Experiment

Factory function to create a [SOMA experiment](SOMAExperiment.md) for
writing (lifecycle: maturing).

## Usage

``` r
SOMAExperimentCreate(
  uri,
  ingest_mode = c("write", "resume"),
  platform_config = NULL,
  tiledbsoma_ctx = NULL,
  tiledb_timestamp = NULL
)
```

## Arguments

- uri:

  URI for the TileDB object.

- ingest_mode:

  Ingestion mode when creating the TileDB object; choose from:

  - “`write`”: create a new TileDB object and error if it already
    exists.

  - “`resume`”: attempt to create a new TileDB object; if it already
    exists, simply open it for writing.

- platform_config:

  Optional platform configuration.

- tiledbsoma_ctx:

  Optional SOMATileDBContext.

- tiledb_timestamp:

  Optional Datetime (POSIXct) for TileDB timestamp.

## Value

A new [SOMA experiment](SOMAExperiment.md) stored at `uri` opened for
writing.

## Examples

``` r
uri <- withr::local_tempfile(pattern = "soma-experiment")
obs <- data.frame(
  soma_joinid = bit64::seq.integer64(0L, 99L),
  obs_id = paste0("cell_", seq_len(100L))
)
sch <- arrow::infer_schema(obs)

(exp <- SOMAExperimentCreate(uri))
#> <SOMAExperiment>
#>   uri: /tmp/Rtmpqg2whH/soma-experiment2b787540efb6
sdf <- exp$add_new_dataframe(
  "obs",
  sch,
  "soma_joinid",
  list(soma_joinid = c(0, 100))
)
sdf$write(arrow::as_arrow_table(obs, schema = sch))
sdf$close()
exp$close()

(exp <- SOMAExperimentOpen(uri))
#> <SOMAExperiment>
#>   uri: /tmp/Rtmpqg2whH/soma-experiment2b787540efb6
exp$obs
#> <SOMADataFrame>
#>   uri: file:///tmp/Rtmpqg2whH/soma-experiment2b787540efb6/obs
#>   dimensions: soma_joinid 
#>   attributes: obs_id 
```
