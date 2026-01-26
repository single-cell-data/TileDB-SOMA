# Open a SOMA Data Frame

Factory function to open a [SOMA data frame](SOMADataFrame.md) for
reading (lifecycle: maturing).

## Usage

``` r
SOMADataFrameOpen(
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
  current time.

- context:

  Optional `SOMAContext` object used for TileDB operations. If a context
  is not provided, then the default context will be used. Call
  `set_default_context` once before other SOMA operations to configure
  the default context.

## Value

A [SOMA data frame](SOMADataFrame.md) stored at `uri` opened in mode
`mode`.

## Examples

``` r
uri <- withr::local_tempfile(pattern = "soma-data-frame")
df <- data.frame(
  soma_joinid = bit64::seq.integer64(0L, 99L),
  group = sample(factor(c("g1", "g2")), size = 100L, replace = TRUE),
  nCount = stats::rbinom(100L, 10L, 0.3)
)
(sch <- arrow::infer_schema(df))
#> Schema
#> soma_joinid: int64
#> group: dictionary<values=string, indices=int8>
#> nCount: int32
(sdf <- SOMADataFrameCreate(uri, sch, domain = list(soma_joinid = c(0, 100))))
#> <SOMADataFrame>
#>   uri: /tmp/RtmpbAgXbM/soma-data-frame284668645d47
#>   dimensions: soma_joinid 
#>   attributes: group, nCount 
sdf$write(arrow::as_arrow_table(df, schema = sch))
sdf$close()

(sdf <- SOMADataFrameOpen(uri))
#> <SOMADataFrame>
#>   uri: /tmp/RtmpbAgXbM/soma-data-frame284668645d47
#>   dimensions: soma_joinid 
#>   attributes: group, nCount 
head(as.data.frame(sdf$read()$concat()))
#>   soma_joinid group nCount
#> 1           0    g2      0
#> 2           1    g1      4
#> 3           2    g2      6
#> 4           3    g2      3
#> 5           4    g1      2
#> 6           5    g2      1
```
