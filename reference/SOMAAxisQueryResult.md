# `SOMAExperiment` Axis Query Result

Access [`SOMAExperimentAxisQuery`](SOMAExperimentAxisQuery.md) results.

## Active bindings

- `obs`:

  [`arrow::Table`](https://arrow.apache.org/docs/r/reference/Table-class.html)
  containing `obs` query slice.

- `var`:

  [`arrow::Table`](https://arrow.apache.org/docs/r/reference/Table-class.html)
  containing `var` query slice. `measurement_name`.

- `X_layers`:

  named list of
  [`arrow::Table`](https://arrow.apache.org/docs/r/reference/Table-class.html)s
  for each `X` layer.

## Methods

### Public methods

- [`SOMAAxisQueryResult$new()`](#method-SOMAAxisQueryResult-new)

- [`SOMAAxisQueryResult$clone()`](#method-SOMAAxisQueryResult-clone)

------------------------------------------------------------------------

### Method `new()`

Create a new `SOMAAxisQueryResult` object.

#### Usage

    SOMAAxisQueryResult$new(obs, var, X_layers)

#### Arguments

- `obs, var`:

  [`arrow::Table`](https://arrow.apache.org/docs/r/reference/Table-class.html)
  containing `obs` or `var` query slice.

- `X_layers`:

  named list of
  [`arrow::Table`](https://arrow.apache.org/docs/r/reference/Table-class.html)s,
  one for each `X` layer.

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    SOMAAxisQueryResult$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
