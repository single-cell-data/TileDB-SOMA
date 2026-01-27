# SOMADataFrame

A SOMA data frame is a multi-column table that must contain a column
called “`soma_joinid`” of type `int64`, which contains a unique value
for each row and is intended to act as a join key for other objects,
such as [`SOMASparseNDArray`](SOMASparseNDArray.md) (lifecycle:
maturing).

## Super classes

[`tiledbsoma::SOMAObject`](SOMAObject.md) -\>
[`tiledbsoma::SOMAArrayBase`](SOMAArrayBase.md) -\> `SOMADataFrame`

## Methods

### Public methods

- [`SOMADataFrame$create()`](#method-SOMADataFrame-create)

- [`SOMADataFrame$write()`](#method-SOMADataFrame-write)

- [`SOMADataFrame$read()`](#method-SOMADataFrame-read)

- [`SOMADataFrame$update()`](#method-SOMADataFrame-update)

- [`SOMADataFrame$levels()`](#method-SOMADataFrame-levels)

- [`SOMADataFrame$shape()`](#method-SOMADataFrame-shape)

- [`SOMADataFrame$maxshape()`](#method-SOMADataFrame-maxshape)

- [`SOMADataFrame$domain()`](#method-SOMADataFrame-domain)

- [`SOMADataFrame$maxdomain()`](#method-SOMADataFrame-maxdomain)

- [`SOMADataFrame$tiledbsoma_has_upgraded_domain()`](#method-SOMADataFrame-tiledbsoma_has_upgraded_domain)

- [`SOMADataFrame$tiledbsoma_resize_soma_joinid_shape()`](#method-SOMADataFrame-tiledbsoma_resize_soma_joinid_shape)

- [`SOMADataFrame$tiledbsoma_upgrade_domain()`](#method-SOMADataFrame-tiledbsoma_upgrade_domain)

- [`SOMADataFrame$change_domain()`](#method-SOMADataFrame-change_domain)

- [`SOMADataFrame$clone()`](#method-SOMADataFrame-clone)

Inherited methods

- [`tiledbsoma::SOMAObject$class()`](SOMAObject.html#method-class)
- [`tiledbsoma::SOMAObject$exists()`](SOMAObject.html#method-exists)
- [`tiledbsoma::SOMAObject$get_metadata()`](SOMAObject.html#method-get_metadata)
- [`tiledbsoma::SOMAObject$initialize()`](SOMAObject.html#method-initialize)
- [`tiledbsoma::SOMAObject$is_open()`](SOMAObject.html#method-is_open)
- [`tiledbsoma::SOMAObject$mode()`](SOMAObject.html#method-mode)
- [`tiledbsoma::SOMAObject$reopen()`](SOMAObject.html#method-reopen)
- [`tiledbsoma::SOMAObject$set_metadata()`](SOMAObject.html#method-set_metadata)
- [`tiledbsoma::SOMAArrayBase$allows_duplicates()`](SOMAArrayBase.html#method-allows_duplicates)
- [`tiledbsoma::SOMAArrayBase$attributes()`](SOMAArrayBase.html#method-attributes)
- [`tiledbsoma::SOMAArrayBase$attrnames()`](SOMAArrayBase.html#method-attrnames)
- [`tiledbsoma::SOMAArrayBase$close()`](SOMAArrayBase.html#method-close)
- [`tiledbsoma::SOMAArrayBase$colnames()`](SOMAArrayBase.html#method-colnames)
- [`tiledbsoma::SOMAArrayBase$dimensions()`](SOMAArrayBase.html#method-dimensions)
- [`tiledbsoma::SOMAArrayBase$dimnames()`](SOMAArrayBase.html#method-dimnames)
- [`tiledbsoma::SOMAArrayBase$index_column_names()`](SOMAArrayBase.html#method-index_column_names)
- [`tiledbsoma::SOMAArrayBase$is_sparse()`](SOMAArrayBase.html#method-is_sparse)
- [`tiledbsoma::SOMAArrayBase$ndim()`](SOMAArrayBase.html#method-ndim)
- [`tiledbsoma::SOMAArrayBase$non_empty_domain()`](SOMAArrayBase.html#method-non_empty_domain)
- [`tiledbsoma::SOMAArrayBase$open()`](SOMAArrayBase.html#method-open)
- [`tiledbsoma::SOMAArrayBase$print()`](SOMAArrayBase.html#method-print)
- [`tiledbsoma::SOMAArrayBase$schema()`](SOMAArrayBase.html#method-schema)

------------------------------------------------------------------------

### Method `create()`

Create a SOMA data frame (lifecycle: maturing).  
  
**Note**: `$create()` is considered internal and should not be called
directly; use factory functions (eg.
[`SOMADataFrameCreate()`](SOMADataFrameCreate.md)) instead.

#### Usage

    SOMADataFrame$create(
      schema,
      index_column_names = c("soma_joinid"),
      domain = NULL,
      platform_config = NULL
    )

#### Arguments

- `schema`:

  An [Arrow
  schema](https://arrow.apache.org/docs/r/reference/schema.html).

- `index_column_names`:

  A vector of column names to use as user-defined index columns. All
  named columns must exist in the schema, and at least one index column
  name is required.

- `domain`:

  An optional list specifying the domain of each index column. Each slot
  in the list must have its name being the name of an index column, and
  its value being be a length-two vector consisting of the minimum and
  maximum values storable in the index column. For example, if there is
  a single int64-valued index column `soma_joinid`, then `domain` might
  be `list(soma_joinid=c(100, 200))` to indicate that values between 100
  and 200, inclusive, can be stored in that column. If provided, this
  sequence must have the same length as `index_column_names`, and the
  index-column domain will be as specified. Omitting or setting the
  domain to `NULL` is deprecated. See also `change_domain` which allows
  you to expand the domain after create.

- `platform_config`:

  A [platform configuration](PlatformConfig.md) object

#### Returns

Returns `self`.

------------------------------------------------------------------------

### Method [`write()`](https://rdrr.io/r/base/write.html)

Write values to the data frame (lifecycle: maturing).

#### Usage

    SOMADataFrame$write(values)

#### Arguments

- `values`:

  An [Arrow
  table](https://arrow.apache.org/docs/r/reference/Table-class.html) or
  [Arrow record
  batch](https://arrow.apache.org/docs/r/reference/RecordBatch-class.html)
  containing all columns, including any index columns. The schema for
  `values` must match the schema for the data frame.

#### Returns

Invisibly returns `self`.

------------------------------------------------------------------------

### Method `read()`

Read a user-defined subset of data, addressed by the data frame indexing
column, and optionally filtered (lifecycle: maturing).

#### Usage

    SOMADataFrame$read(
      coords = NULL,
      column_names = NULL,
      value_filter = NULL,
      result_order = "auto",
      log_level = "auto"
    )

#### Arguments

- `coords`:

  Optional named list of indices specifying the rows to read; each
  (named) list element corresponds to a dimension of the same name.

- `column_names`:

  Optional character vector of column names to return.

- `value_filter`:

  Optional string containing a logical expression that is used to filter
  the returned values. See
  [`tiledb::parse_query_condition()`](https://tiledb-inc.github.io/TileDB-R/reference/parse_query_condition.html)
  for more information.

- `result_order`:

  Optional order of read results. This can be one of either
  `"ROW_MAJOR, `"COL_MAJOR"`, or `"auto"\` (default).

- `log_level`:

  Optional logging level with default value of “`warn`”.

#### Returns

An [Arrow
table](https://arrow.apache.org/docs/r/reference/Table-class.html) or
[`TableReadIter`](TableReadIter.md)

------------------------------------------------------------------------

### Method [`update()`](https://rdrr.io/r/stats/update.html)

Update (lifecycle: maturing).

#### Usage

    SOMADataFrame$update(values, row_index_name = NULL)

#### Arguments

- `values`:

  A data frame, [Arrow
  table](https://arrow.apache.org/docs/r/reference/Table-class.html), or
  [Arrow record
  batch](https://arrow.apache.org/docs/r/reference/RecordBatch-class.html).

- `row_index_name`:

  An optional scalar character. If provided, and if the `values`
  argument is a data frame with row names, then the row names will be
  extracted and added as a new column to the data frame prior to
  performing the update. The name of this new column will be set to the
  value specified by `row_index_name`.

#### Details

Update the existing `SOMADataFrame` to add or remove columns based on
the input:

- columns present in the current the `SOMADataFrame` but absent from the
  new `values` will be dropped.

- columns absent in current `SOMADataFrame` but present in the new
  `values` will be added.

- any columns present in both will be left alone, with the exception
  that if `values` has a different type for the column, the entire
  update will fail because attribute types cannot be changed.

Furthermore, `values` must contain the same number of rows as the
current `SOMADataFrame`.

#### Returns

Invisibly returns `NULL`

------------------------------------------------------------------------

### Method [`levels()`](https://rdrr.io/r/base/levels.html)

Get the levels for an enumerated (`factor`) column.

#### Usage

    SOMADataFrame$levels(column_names = NULL, simplify = TRUE)

#### Arguments

- `column_names`:

  Optional character vector of column names to pull enumeration levels
  for; defaults to all enumerated columns.

- `simplify`:

  Simplify the result down to a vector or matrix.

#### Returns

If `simplify` returns one of the following:

- a vector of there is only one enumerated column.

- a matrix if there are multiple enumerated columns with the same number
  of levels.

- a named list if there are multiple enumerated columns with differing
  numbers of levels.

Otherwise, returns a named list.

------------------------------------------------------------------------

### Method `shape()`

Retrieve the shape; as `SOMADataFrames` are shapeless, simply raises an
error.

#### Usage

    SOMADataFrame$shape()

#### Returns

None, instead a
[`.NotYetImplemented()`](https://rdrr.io/r/base/notyet.html) error is
raised.

------------------------------------------------------------------------

### Method `maxshape()`

Retrieve the max shape; as `SOMADataFrames` are shapeless, simply raises
an error.

#### Usage

    SOMADataFrame$maxshape()

#### Returns

None, instead a
[`.NotYetImplemented()`](https://rdrr.io/r/base/notyet.html) error is
raised.

------------------------------------------------------------------------

### Method `domain()`

Returns a named list of minimum/maximum pairs, one per index column,
currently storable on each index column of the data frame. These can be
resized up to `maxdomain` (lifecycle: maturing).

#### Usage

    SOMADataFrame$domain()

#### Returns

Named list of minimum/maximum values.

------------------------------------------------------------------------

### Method `maxdomain()`

Returns a named list of minimum/maximum pairs, one per index column,
which are the limits up to which the data frame can have its domain
resized (lifecycle: maturing).

#### Usage

    SOMADataFrame$maxdomain()

#### Returns

Named list of minimum/maximum values.

------------------------------------------------------------------------

### Method `tiledbsoma_has_upgraded_domain()`

Test if the array has the upgraded resizeable domain feature from
TileDB-SOMA 1.15, the array was created with this support, or it has had
`$upgrade_domain()` applied to it (lifecycle: maturing).

#### Usage

    SOMADataFrame$tiledbsoma_has_upgraded_domain()

#### Returns

Returns `TRUE` if the array has the upgraded resizable domain feature;
otherwise, returns `FALSE`.

------------------------------------------------------------------------

### Method `tiledbsoma_resize_soma_joinid_shape()`

Increases the shape of the data frame on the `soma_joinid` index column,
if it indeed is an index column, leaving all other index columns as-is.
If the `soma_joinid` is not an index column, no change is made. This is
a special case of `upgrade_domain()`, but simpler to keystroke, and
handles the most common case for data frame domain expansion. Raises an
error if the data frame doesn't already have a domain; in that case
please call `$tiledbsoma_upgrade_domain()`.

#### Usage

    SOMADataFrame$tiledbsoma_resize_soma_joinid_shape(new_shape)

#### Arguments

- `new_shape`:

  An integer, greater than or equal to 1 + the `soma_joinid` domain
  slot.

#### Returns

Invisibly returns `NULL`

------------------------------------------------------------------------

### Method `tiledbsoma_upgrade_domain()`

Allows you to set the domain of a `SOMADataFrame`, when the
`SOMADataFrame` does not have a domain set yet. The argument must be a
list of pairs of low/high values for the desired domain, one pair per
index column. For string index columns, you must offer the low/high pair
as `c("", "")`, or as `NULL`. If `check_only` is `True`, returns whether
the operation would succeed if attempted, or a reason why it would not.
The domain being requested must be contained within what `$maxdomain()`
returns.

#### Usage

    SOMADataFrame$tiledbsoma_upgrade_domain(new_domain, check_only = FALSE)

#### Arguments

- `new_domain`:

  A named list, keyed by index-column name, with values being
  two-element vectors containing the desired lower and upper bounds for
  the domain.

- `check_only`:

  If true, does not apply the operation, but only reports whether it
  would have succeeded.

#### Returns

If `check_only`, returns the empty string if no error is detected, else
a description of the error. Otherwise, invisibly returns `NULL`

------------------------------------------------------------------------

### Method `change_domain()`

Allows you to set the domain of a `SOMADataFrame`, when the
`SOMADataFrame` already has a domain set yet. The argument must be a
list of pairs of low/high values for the desired domain, one pair per
index column. For string index columns, you must offer the low/high pair
as `c("", "")`, or as `NULL`. If `check_only` is `True`, returns whether
the operation would succeed if attempted, or a reason why it would not.
The return value from `domain` must be contained within the requested
`new_domain`, and the requested `new_domain` must be contained within
the return value from `$maxdomain()` (lifecycle: maturing).

#### Usage

    SOMADataFrame$change_domain(new_domain, check_only = FALSE)

#### Arguments

- `new_domain`:

  A named list, keyed by index-column name, with values being
  two-element vectors containing the desired lower and upper bounds for
  the domain.

- `check_only`:

  If true, does not apply the operation, but only reports whether it
  would have succeeded.

#### Returns

If `check_only`, returns the empty string if no error is detected, else
a description of the error. Otherwise, invisibly returns `NULL`

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    SOMADataFrame$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

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
#>   uri: /tmp/Rtmpfr8mYm/soma-data-frame2a7e36156b76
#>   dimensions: soma_joinid 
#>   attributes: group, nCount 
sdf$write(arrow::as_arrow_table(df, schema = sch))
sdf$close()

(sdf <- SOMADataFrameOpen(uri))
#> <SOMADataFrame>
#>   uri: /tmp/Rtmpfr8mYm/soma-data-frame2a7e36156b76
#>   dimensions: soma_joinid 
#>   attributes: group, nCount 
head(as.data.frame(sdf$read()$concat()))
#>   soma_joinid group nCount
#> 1           0    g2      3
#> 2           1    g1      4
#> 3           2    g2      1
#> 4           3    g2      2
#> 5           4    g1      2
#> 6           5    g1      4
```
