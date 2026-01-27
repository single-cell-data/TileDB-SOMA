# SOMA Blockwise Read Iterator for Arrow Tables

Class that allows for blockwise read iteration of SOMA reads as Arrow
[`Table`](https://arrow.apache.org/docs/r/reference/Table-class.html)`s`.

## Super classes

[`tiledbsoma::ReadIter`](ReadIter.md) -\>
[`tiledbsoma::BlockwiseReadIterBase`](BlockwiseReadIterBase.md) -\>
`BlockwiseTableReadIter`

## Methods

### Public methods

- [`BlockwiseTableReadIter$concat()`](#method-BlockwiseTableReadIter-concat)

- [`BlockwiseTableReadIter$clone()`](#method-BlockwiseTableReadIter-clone)

Inherited methods

- [`tiledbsoma::BlockwiseReadIterBase$initialize()`](BlockwiseReadIterBase.html#method-initialize)
- [`tiledbsoma::BlockwiseReadIterBase$read_complete()`](BlockwiseReadIterBase.html#method-read_complete)
- [`tiledbsoma::BlockwiseReadIterBase$read_next()`](BlockwiseReadIterBase.html#method-read_next)

------------------------------------------------------------------------

### Method `concat()`

Concatenate the remainder of the blockwise iterator.

#### Usage

    BlockwiseTableReadIter$concat()

#### Returns

An Arrow Table with the remainder of the iterator.

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    BlockwiseTableReadIter$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
dir <- withr::local_tempfile(pattern = "blockwise-table")
dir.create(dir, recursive = TRUE)
(exp <- load_dataset("soma-exp-pbmc-small", dir))
#> <SOMAExperiment>
#>   uri: /tmp/Rtmpfr8mYm/blockwise-table2a7e51c7ba34/soma-exp-pbmc-small
qry <- exp$axis_query("RNA")
xqry <- qry$X("data")

iter <- xqry$blockwise(0L, size = 20L, reindex_disable_on_axis = TRUE)$tables()
stopifnot(inherits(iter, "BlockwiseTableReadIter"))

while (!iter$read_complete()) {
  block <- iter$read_next()
}
```
