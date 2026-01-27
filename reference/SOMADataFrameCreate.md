# Create a SOMA Data Frame

Factory function to create a [SOMA data frame](SOMADataFrame.md) for
writing (lifecycle: maturing).

## Usage

``` r
SOMADataFrameCreate(
  uri,
  schema,
  index_column_names = c("soma_joinid"),
  domain = NULL,
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

- schema:

  Arrow schema argument for the [SOMA dataframe](SOMADataFrame.md).

- index_column_names:

  A vector of column names to use as user-defined index columns; all
  named columns must exist in the schema, and at least one index column
  name is required.

- domain:

  An optional list of 2-element vectors specifying the domain of each
  index column. Each vector should be a pair consisting of the minimum
  and maximum values storable in the index column. For example, if there
  is a single int64-valued index column, then `domain` might be
  `c(100, 200)` to indicate that values between 100 and 200, inclusive,
  can be stored in that column. If provided, this list must have the
  same length as `index_column_names`, and the index-column domain will
  be as specified. If omitted entirely, or if `NULL` in a given
  dimension, the corresponding index-column domain will use the minimum
  and maximum possible values for the column's datatype. This makes a
  SOMA data frame growable.

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

A new [SOMA data frame](SOMADataFrame.md) stored at `uri` opened for
writing.

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
#>   uri: /tmp/Rtmpfr8mYm/soma-data-frame2a7e5df9cde5
#>   dimensions: soma_joinid 
#>   attributes: group, nCount 
sdf$write(arrow::as_arrow_table(df, schema = sch))
sdf$close()

(sdf <- SOMADataFrameOpen(uri))
#> <SOMADataFrame>
#>   uri: /tmp/Rtmpfr8mYm/soma-data-frame2a7e5df9cde5
#>   dimensions: soma_joinid 
#>   attributes: group, nCount 
head(as.data.frame(sdf$read()$concat()))
#>   soma_joinid group nCount
#> 1           0    g1      1
#> 2           1    g1      4
#> 3           2    g1      2
#> 4           3    g2      0
#> 5           4    g1      1
#> 6           5    g1      4
```
