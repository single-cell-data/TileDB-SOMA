# SOMA Read Iterator Over Arrow Tables

`TableReadIter` is a class that allows for iteration over a reads on
[SOMASparseNDArray](SOMASparseNDArray.md) and
[SOMADataFrame](SOMADataFrame.md). Iteration chunks are retrieved as an
Arrow
[Table](https://arrow.apache.org/docs/r/reference/Table-class.html).

## Super class

[`tiledbsoma::ReadIter`](ReadIter.md) -\> `TableReadIter`

## Methods

### Public methods

- [`TableReadIter$concat()`](#method-TableReadIter-concat)

- [`TableReadIter$clone()`](#method-TableReadIter-clone)

Inherited methods

- [`tiledbsoma::ReadIter$initialize()`](ReadIter.html#method-initialize)
- [`tiledbsoma::ReadIter$read_complete()`](ReadIter.html#method-read_complete)
- [`tiledbsoma::ReadIter$read_next()`](ReadIter.html#method-read_next)

------------------------------------------------------------------------

### Method `concat()`

Concatenate remainder of iterator.

#### Usage

    TableReadIter$concat()

#### Returns

An Arrow
[`Table`](https://arrow.apache.org/docs/r/reference/Table-class.html).

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    TableReadIter$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
dir <- withr::local_tempfile(pattern = "table-iter")
dir.create(dir, recursive = TRUE)
(exp <- load_dataset("soma-exp-pbmc-small", dir))
#> <SOMAExperiment>
#>   uri: /tmp/Rtmpfr8mYm/table-iter2a7e7f1b400e/soma-exp-pbmc-small
qry <- exp$axis_query("RNA")
xqry <- qry$X("data")

iter <- xqry$tables()
stopifnot(inherits(iter, "TableReadIter"))

while (!iter$read_complete()) {
  block <- iter$read_next()
}
```
