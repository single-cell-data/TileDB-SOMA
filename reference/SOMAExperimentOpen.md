# Open SOMA Experiment

Factory function to open a [SOMA experiment](SOMAExperiment.md) for
reading (lifecycle: maturing).

## Usage

``` r
SOMAExperimentOpen(
  uri,
  mode = "READ",
  platform_config = NULL,
  tiledbsoma_ctx = NULL,
  tiledb_timestamp = NULL,
  context = NULL
)
```

## Arguments

- uri:

  URI for the TileDB object.

- mode:

  One of “`READ`” or “`WRITE`”.

- platform_config:

  Optional platform configuration.

- tiledbsoma_ctx:

  Optional (DEPRECATED) SOMATileDBContext.

- tiledb_timestamp:

  Optional Datetime (POSIXct) for TileDB timestamp; defaults to the
  current time. If not `NULL`, all members accessed through the
  collection inherit the timestamp.

- context:

  Optional `SOMAContext` object used for TileDB operations. If a context
  is not provided, then the default context will be used. Call
  `set_default_context` once before other SOMA operations to configure
  the default context.

## Value

A [SOMA experiment](SOMAExperiment.md) stored at `uri` opened in mode
`mode`.

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
#>   uri: /tmp/Rtmpfr8mYm/soma-experiment2a7e515d850d
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
#>   uri: /tmp/Rtmpfr8mYm/soma-experiment2a7e515d850d
exp$obs
#> <SOMADataFrame>
#>   uri: file:///tmp/Rtmpfr8mYm/soma-experiment2a7e515d850d/obs
#>   dimensions: soma_joinid 
#>   attributes: obs_id 
```
