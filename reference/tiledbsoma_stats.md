# TileDB SOMA Statistics

These functions expose the TileDB Core functionality for performance
measurements and statistics

## Usage

``` r
tiledbsoma_stats_enable()

tiledbsoma_stats_disable()

tiledbsoma_stats_reset()

tiledbsoma_stats_dump()

tiledbsoma_stats_show()
```

## Value

`tiledbsoma_stats_show()`: a single-length character vector with the
TileDB statistics encoded in JSON format

All other functions invisibly return `NULL`

## Details

- `tiledbsoma_stats_enable()`/`tiledbsoma_stats_disable()`: Enable and
  disable TielDB's internal statistics

- `tiledbsoma_stats_reset()`: Reset all statistics to `0`

- `tiledbsoma_stats_dump()`: Dump all statistics as a JSON string

- `tiledbsoma_stats_show()`: Pretty-print the JSON statistics
