# `SOMAExperiment` Axis Query Result Indexer

Index `obs`/`var` soma_joinids for a given query result.

Retrieve the index of the given `obs` or `var` coordinates in the query
result. Coordinates outside of the query result will return
[`arrow::null()`](https://arrow.apache.org/docs/r/reference/data-type.html).

## Methods

### Public methods

- [`SOMAAxisIndexer$new()`](#method-SOMAAxisIndexer-new)

- [`SOMAAxisIndexer$by_obs()`](#method-SOMAAxisIndexer-by_obs)

- [`SOMAAxisIndexer$by_var()`](#method-SOMAAxisIndexer-by_var)

- [`SOMAAxisIndexer$clone()`](#method-SOMAAxisIndexer-clone)

------------------------------------------------------------------------

### Method `new()`

Create a new `SOMAAxisIndexer` object.

#### Usage

    SOMAAxisIndexer$new(query)

#### Arguments

- `query`:

  The [`SOMAExperimentAxisQuery`](SOMAExperimentAxisQuery.md) object to
  build indices for.

------------------------------------------------------------------------

### Method `by_obs()`

Get the index of the given `obs` coordinates.

#### Usage

    SOMAAxisIndexer$by_obs(coords)

#### Arguments

- `coords`:

  vector or
  [`arrow::Array`](https://arrow.apache.org/docs/r/reference/array-class.html)
  of numeric coordinates.

------------------------------------------------------------------------

### Method `by_var()`

Get the index of the given `var` coordinates.

#### Usage

    SOMAAxisIndexer$by_var(coords)

#### Arguments

- `coords`:

  vector or
  [`arrow::Array`](https://arrow.apache.org/docs/r/reference/array-class.html)
  of numeric coordinates.

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    SOMAAxisIndexer$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
