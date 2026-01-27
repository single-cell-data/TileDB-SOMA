# Create a SOMA Measurement

Factory function to create a [SOMA measurement](SOMAMeasurement.md) for
writing (lifecycle: maturing).

## Usage

``` r
SOMAMeasurementCreate(
  uri,
  ingest_mode = c("write", "resume"),
  platform_config = NULL,
  tiledbsoma_ctx = NULL,
  tiledb_timestamp = NULL,
  context = NULL
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

  Optional (DEPRECATED) SOMATileDBContext.

- tiledb_timestamp:

  Optional Datetime (POSIXct) for TileDB timestamp.

- context:

  Optional `SOMAContext` object used for TileDB operations. If a context
  is not provided, then the default context will be used. Call
  `set_default_context` once before other SOMA operations to configure
  the default context.

## Value

A new [SOMA measurement](SOMAMeasurement.md) stored at `uri` opened for
writing.

## Examples

``` r
uri <- withr::local_tempfile(pattern = "soma-measurement")
var <- data.frame(
  soma_joinid = bit64::seq.integer64(0L, 99L),
  var_id = paste0("feature_", seq_len(100L))
)
sch <- arrow::infer_schema(var)

(ms <- SOMAMeasurementCreate(uri))
#> <SOMAMeasurement>
#>   uri: /tmp/Rtmpfr8mYm/soma-measurement2a7ec1bed4a
sdf <- ms$add_new_dataframe(
  "var",
  sch,
  "soma_joinid",
  list(soma_joinid = c(0, 100))
)
sdf$write(arrow::as_arrow_table(var, schema = sch))
sdf$close()
ms$close()

(ms <- SOMAMeasurementOpen(uri))
#> <SOMAMeasurement>
#>   uri: /tmp/Rtmpfr8mYm/soma-measurement2a7ec1bed4a
ms$var
#> <SOMADataFrame>
#>   uri: file:///tmp/Rtmpfr8mYm/soma-measurement2a7ec1bed4a/var
#>   dimensions: soma_joinid 
#>   attributes: var_id 
```
