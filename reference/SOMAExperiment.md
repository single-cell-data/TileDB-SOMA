# SOMA Experiment

`SOMAExperiment` is a specialized [`SOMACollection`](SOMACollection.md),
representing one or more modes of measurement across a single collection
of cells (aka a “multimodal dataset”) with pre-defined fields: `obs` and
`ms` (see *Active Bindings* below for details) (lifecycle: maturing).

## Adding new objects to a collection

The `SOMAExperiment` class provides a number of type-specific methods
for adding new a object to the collection, such as
`add_new_sparse_ndarray()` and `add_new_dataframe()`. These methods will
create the new object and add it as member of the `SOMAExperiment`. The
new object will always inherit the parent context (see
[`SOMATileDBContext`](SOMATileDBContext.md)) and, by default, its
platform configuration (see [`PlatformConfig`](PlatformConfig.md)).
However, the user can override the default platform configuration by
passing a custom configuration to the `platform_config` argument.

## Carrara (TileDB v3) behavior

When working with Carrara URIs (`tiledb://workspace/teamspace/...`),
child objects created at a URI nested under a parent collection are
**automatically added** as members of the parent. This means:

- You do not need to call `add_new_collection()` after creating a child
  at a nested URI—the child is already a member.

- For backward compatibility, calling `add_new_collection()` on an
  already-registered child is a **no-op** and will not cause an error.

- The member name must match the relative URI segment (e.g., creating at
  `parent_uri/child` automatically adds the child with key `"child"`).

## Super classes

[`tiledbsoma::SOMAObject`](SOMAObject.md) -\>
[`tiledbsoma::SOMACollectionBase`](SOMACollectionBase.md) -\>
`SOMAExperiment`

## Active bindings

- `obs`:

  A [`SOMADataFrame`](SOMADataFrame.md) containing primary annotations
  on the observation axis. The contents of the `soma_joinid` column
  define the observation index domain, `obs_id`. All observations for
  the `SOMAExperiment` must be defined in this data frame.

- `ms`:

  A [`SOMACollection`](SOMACollection.md) of named
  [`SOMAMeasurement`](SOMAMeasurement.md)`s`.

## Methods

### Public methods

- [`SOMAExperiment$axis_query()`](#method-SOMAExperiment-axis_query)

- [`SOMAExperiment$update_obs()`](#method-SOMAExperiment-update_obs)

- [`SOMAExperiment$update_var()`](#method-SOMAExperiment-update_var)

- [`SOMAExperiment$clone()`](#method-SOMAExperiment-clone)

Inherited methods

- [`tiledbsoma::SOMAObject$class()`](SOMAObject.html#method-class)
- [`tiledbsoma::SOMAObject$exists()`](SOMAObject.html#method-exists)
- [`tiledbsoma::SOMAObject$get_metadata()`](SOMAObject.html#method-get_metadata)
- [`tiledbsoma::SOMAObject$initialize()`](SOMAObject.html#method-initialize)
- [`tiledbsoma::SOMAObject$is_open()`](SOMAObject.html#method-is_open)
- [`tiledbsoma::SOMAObject$mode()`](SOMAObject.html#method-mode)
- [`tiledbsoma::SOMAObject$reopen()`](SOMAObject.html#method-reopen)
- [`tiledbsoma::SOMACollectionBase$add_new_collection()`](SOMACollectionBase.html#method-add_new_collection)
- [`tiledbsoma::SOMACollectionBase$add_new_dataframe()`](SOMACollectionBase.html#method-add_new_dataframe)
- [`tiledbsoma::SOMACollectionBase$add_new_dense_ndarray()`](SOMACollectionBase.html#method-add_new_dense_ndarray)
- [`tiledbsoma::SOMACollectionBase$add_new_sparse_ndarray()`](SOMACollectionBase.html#method-add_new_sparse_ndarray)
- [`tiledbsoma::SOMACollectionBase$close()`](SOMACollectionBase.html#method-close)
- [`tiledbsoma::SOMACollectionBase$create()`](SOMACollectionBase.html#method-create)
- [`tiledbsoma::SOMACollectionBase$get()`](SOMACollectionBase.html#method-get)
- [`tiledbsoma::SOMACollectionBase$length()`](SOMACollectionBase.html#method-length)
- [`tiledbsoma::SOMACollectionBase$names()`](SOMACollectionBase.html#method-names)
- [`tiledbsoma::SOMACollectionBase$open()`](SOMACollectionBase.html#method-open)
- [`tiledbsoma::SOMACollectionBase$print()`](SOMACollectionBase.html#method-print)
- [`tiledbsoma::SOMACollectionBase$remove()`](SOMACollectionBase.html#method-remove)
- [`tiledbsoma::SOMACollectionBase$set()`](SOMACollectionBase.html#method-set)
- [`tiledbsoma::SOMACollectionBase$set_metadata()`](SOMACollectionBase.html#method-set_metadata)

------------------------------------------------------------------------

### Method `axis_query()`

Subset and extract data from a [`SOMAMeasurement`](SOMAMeasurement.md)
by querying the `obs`/`var` axes.

#### Usage

    SOMAExperiment$axis_query(measurement_name, obs_query = NULL, var_query = NULL)

#### Arguments

- `measurement_name`:

  The name of the measurement to query.

- `obs_query, var_query`:

  An [`SOMAAxisQuery`](SOMAAxisQuery.md) object for the obs/var axis.

#### Returns

A [`SOMAExperimentAxisQuery`](SOMAExperimentAxisQuery.md) object.

------------------------------------------------------------------------

### Method `update_obs()`

Update the `obs` data frame to add or remove columns. See
[`SOMADataFrame$update()`](SOMADataFrame.md) for more details.

#### Usage

    SOMAExperiment$update_obs(values, row_index_name = NULL)

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

------------------------------------------------------------------------

### Method `update_var()`

Update the `var` data frame to add or remove columns. See
[`SOMADataFrame$update()`](SOMADataFrame.md) for more details.

#### Usage

    SOMAExperiment$update_var(values, measurement_name, row_index_name = NULL)

#### Arguments

- `values`:

  A data frame, [Arrow
  table](https://arrow.apache.org/docs/r/reference/Table-class.html), or
  [Arrow record
  batch](https://arrow.apache.org/docs/r/reference/RecordBatch-class.html).

- `measurement_name`:

  The name of the [`SOMAMeasurement`](SOMAMeasurement.md) whose `var`
  will be updated.

- `row_index_name`:

  An optional scalar character. If provided, and if the `values`
  argument is a data frame with row names, then the row names will be
  extracted and added as a new column to the data frame prior to
  performing the update. The name of this new column will be set to the
  value specified by `row_index_name`.

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    SOMAExperiment$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
uri <- withr::local_tempfile(pattern = "soma-experiment")
obs <- data.frame(
  soma_joinid = bit64::seq.integer64(0L, 99L),
  obs_id = paste0("cell_", seq_len(100L))
)
sch <- arrow::infer_schema(obs)

(exp <- SOMAExperimentCreate(uri))
#> <SOMAExperiment>
#>   uri: /tmp/Rtmpfr8mYm/soma-experiment2a7e4e423aa1
sdf <- exp$add_new_dataframe(
  "obs",
  sch,
  "soma_joinid",
  list(soma_joinid = c(0, 100))
)
sdf$write(arrow::as_arrow_table(obs, schema = sch))
sdf$close()
exp$close()

(exp <- SOMAExperimentOpen(uri))
#> <SOMAExperiment>
#>   uri: /tmp/Rtmpfr8mYm/soma-experiment2a7e4e423aa1
exp$obs
#> <SOMADataFrame>
#>   uri: file:///tmp/Rtmpfr8mYm/soma-experiment2a7e4e423aa1/obs
#>   dimensions: soma_joinid 
#>   attributes: obs_id 
```
