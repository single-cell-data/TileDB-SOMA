# SOMA Axis Query

Construct a single-axis query object with a combination of coordinates
and/or value filters for use with
[`SOMAExperimentAxisQuery`](SOMAExperimentAxisQuery.md). (lifecycle:
maturing)

Per dimension, the `SOMAAxisQuery` can have value of:

- None (i.e., `coords = NULL` and `value_filter = NULL`) - read all
  values

- Coordinates - a set of coordinates on the axis dataframe index,
  expressed in any type or format supported by
  [`SOMADataFrame`](SOMADataFrame.md)'s `read()` method.

- A SOMA `value_filter` across columns in the axis dataframe, expressed
  as string

- Or, a combination of coordinates and value filter.

## See also

[`tiledb::parse_query_condition()`](https://tiledb-inc.github.io/TileDB-R/reference/parse_query_condition.html)
for more information about valid value filters.

## Public fields

- `coords`:

  The coordinates for the query.

- `value_filter`:

  The value filter for the query.

## Methods

### Public methods

- [`SOMAAxisQuery$new()`](#method-SOMAAxisQuery-new)

- [`SOMAAxisQuery$clone()`](#method-SOMAAxisQuery-clone)

------------------------------------------------------------------------

### Method `new()`

Create a new `SOMAAxisQuery` object.

#### Usage

    SOMAAxisQuery$new(value_filter = NULL, coords = NULL)

#### Arguments

- `value_filter`:

  Optional string containing a logical expression that is used to filter
  the returned values.

- `coords`:

  Optional indices specifying the rows to read: either a vector of the
  appropriate type or a named list of vectors corresponding to each
  dimension.

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    SOMAAxisQuery$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
