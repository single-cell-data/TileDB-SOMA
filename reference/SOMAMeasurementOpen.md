# Open SOMA Measurement

Factory function to open a [SOMA measurement](SOMAMeasurement.md) for
reading (lifecycle: maturing).

## Usage

``` r
SOMAMeasurementOpen(
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

A [SOMA measurement](SOMAMeasurement.md) stored at `uri` opened in mode
`mode`.

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
#>   uri: /tmp/RtmpbAgXbM/soma-measurement2846363a5fe4
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
#>   uri: /tmp/RtmpbAgXbM/soma-measurement2846363a5fe4
ms$var
#> <SOMADataFrame>
#>   uri: file:///tmp/RtmpbAgXbM/soma-measurement2846363a5fe4/var
#>   dimensions: soma_joinid 
#>   attributes: var_id 
```
