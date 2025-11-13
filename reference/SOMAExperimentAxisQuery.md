# `SOMAExperiment` Axis Query

Perform an axis-based query against a
[`SOMAExperiment`](SOMAExperiment.md).

`SOMAExperimentAxisQuery` allows easy selection and extraction of data
from a single [`SOMAMeasurement`](SOMAMeasurement.md) in a
[`SOMAExperiment`](SOMAExperiment.md), by `obs`/`var` (axis) coordinates
and/or value filter. The primary use for this class is slicing
[`SOMAExperiment`](SOMAExperiment.md) `X` layers by `obs` or `var` value
and/or coordinates. (lifecycle: maturing)

### X Layer Support

Slicing on [`SOMASparseNDArray`](SOMASparseNDArray.md) `X` matrices is
supported; slicing on [`SOMADenseNDArray`](SOMADenseNDArray.md) is not
supported at this time.

### Result Size

`SOMAExperimentAxisQuery` query class assumes it can store the full
result of both axis dataframe queries in memory, and only provides
incremental access to the underlying X NDArray. Accessors such as
`n_obs` and `n_vars` codify this in the class.

## Active bindings

- `experiment`:

  The parent [`SOMAExperiment`](SOMAExperiment.md) object.

- `indexer`:

  The [`SOMAAxisIndexer`](SOMAAxisIndexer.md) object.

- `obs_query`:

  The `obs` [`SOMAAxisQuery`](SOMAAxisQuery.md) object.

- `var_query`:

  The `var` [`SOMAAxisQuery`](SOMAAxisQuery.md) object.

- `n_obs`:

  The number of `obs` axis query results.

- `n_vars`:

  The number of `var` axis query results.

- `obs_df`:

  The `obs` [`SOMADataFrame`](SOMADataFrame.md) object.

- `var_df`:

  The `var` [`SOMADataFrame`](SOMADataFrame.md) object for the specified
  `measurement_name`.

- `ms`:

  The [`SOMAMeasurement`](SOMAMeasurement.md) object for the specified
  `measurement_name`.

## Methods

### Public methods

- [`SOMAExperimentAxisQuery$new()`](#method-SOMAExperimentAxisQuery-new)

- [`SOMAExperimentAxisQuery$obs()`](#method-SOMAExperimentAxisQuery-obs)

- [`SOMAExperimentAxisQuery$var()`](#method-SOMAExperimentAxisQuery-var)

- [`SOMAExperimentAxisQuery$obs_joinids()`](#method-SOMAExperimentAxisQuery-obs_joinids)

- [`SOMAExperimentAxisQuery$var_joinids()`](#method-SOMAExperimentAxisQuery-var_joinids)

- [`SOMAExperimentAxisQuery$X()`](#method-SOMAExperimentAxisQuery-X)

- [`SOMAExperimentAxisQuery$obsm()`](#method-SOMAExperimentAxisQuery-obsm)

- [`SOMAExperimentAxisQuery$varm()`](#method-SOMAExperimentAxisQuery-varm)

- [`SOMAExperimentAxisQuery$obsp()`](#method-SOMAExperimentAxisQuery-obsp)

- [`SOMAExperimentAxisQuery$varp()`](#method-SOMAExperimentAxisQuery-varp)

- [`SOMAExperimentAxisQuery$read()`](#method-SOMAExperimentAxisQuery-read)

- [`SOMAExperimentAxisQuery$to_sparse_matrix()`](#method-SOMAExperimentAxisQuery-to_sparse_matrix)

- [`SOMAExperimentAxisQuery$to_seurat()`](#method-SOMAExperimentAxisQuery-to_seurat)

- [`SOMAExperimentAxisQuery$to_seurat_assay()`](#method-SOMAExperimentAxisQuery-to_seurat_assay)

- [`SOMAExperimentAxisQuery$to_seurat_reduction()`](#method-SOMAExperimentAxisQuery-to_seurat_reduction)

- [`SOMAExperimentAxisQuery$to_seurat_graph()`](#method-SOMAExperimentAxisQuery-to_seurat_graph)

- [`SOMAExperimentAxisQuery$to_single_cell_experiment()`](#method-SOMAExperimentAxisQuery-to_single_cell_experiment)

- [`SOMAExperimentAxisQuery$clone()`](#method-SOMAExperimentAxisQuery-clone)

------------------------------------------------------------------------

### Method `new()`

Create a new `SOMAExperimentAxisQuery` object.

#### Usage

    SOMAExperimentAxisQuery$new(
      experiment,
      measurement_name,
      obs_query = NULL,
      var_query = NULL
    )

#### Arguments

- `experiment`:

  A [`SOMAExperiment`](SOMAExperiment.md) object.

- `measurement_name`:

  The name of the measurement to query.

- `obs_query, var_query`:

  An [`SOMAAxisQuery`](SOMAAxisQuery.md) object for the obs/var axis.

------------------------------------------------------------------------

### Method `obs()`

Retrieve obs [TableReadIter](TableReadIter.md)

#### Usage

    SOMAExperimentAxisQuery$obs(column_names = NULL)

#### Arguments

- `column_names`:

  A character vector of column names to retrieve

------------------------------------------------------------------------

### Method [`var()`](https://rdrr.io/r/stats/cor.html)

Retrieve var
[`arrow::Table`](https://arrow.apache.org/docs/r/reference/Table-class.html)

#### Usage

    SOMAExperimentAxisQuery$var(column_names = NULL)

#### Arguments

- `column_names`:

  A character vector of column names to retrieve

------------------------------------------------------------------------

### Method `obs_joinids()`

Retrieve `soma_joinids` as an
[`arrow::Array`](https://arrow.apache.org/docs/r/reference/array-class.html)
for `obs`.

#### Usage

    SOMAExperimentAxisQuery$obs_joinids()

------------------------------------------------------------------------

### Method `var_joinids()`

Retrieve `soma_joinids` as an
[`arrow::Array`](https://arrow.apache.org/docs/r/reference/array-class.html)
for `var`.

#### Usage

    SOMAExperimentAxisQuery$var_joinids()

------------------------------------------------------------------------

### Method `X()`

Retrieves an `X` layer as a
[SOMASparseNDArrayRead](SOMASparseNDArrayRead.md)

#### Usage

    SOMAExperimentAxisQuery$X(layer_name)

#### Arguments

- `layer_name`:

  The name of the layer to retrieve.

------------------------------------------------------------------------

### Method `obsm()`

Retrieves an `obsm` layer as a
[`SOMASparseNDArrayRead`](SOMASparseNDArrayRead.md)

#### Usage

    SOMAExperimentAxisQuery$obsm(layer_name)

#### Arguments

- `layer_name`:

  The name of the layer to retrieve

------------------------------------------------------------------------

### Method `varm()`

Retrieves a `varm` layer as a
[`SOMASparseNDArrayRead`](SOMASparseNDArrayRead.md)

#### Usage

    SOMAExperimentAxisQuery$varm(layer_name)

#### Arguments

- `layer_name`:

  The name of the layer to retrieve

------------------------------------------------------------------------

### Method `obsp()`

Retrieves an `obsp` layer as a
[`SOMASparseNDArrayRead`](SOMASparseNDArrayRead.md)

#### Usage

    SOMAExperimentAxisQuery$obsp(layer_name)

#### Arguments

- `layer_name`:

  The name of the layer to retrieve

------------------------------------------------------------------------

### Method `varp()`

Retrieves a `varp` layer as a
[`SOMASparseNDArrayRead`](SOMASparseNDArrayRead.md)

#### Usage

    SOMAExperimentAxisQuery$varp(layer_name)

#### Arguments

- `layer_name`:

  The name of the layer to retrieve

------------------------------------------------------------------------

### Method `read()`

Reads the entire query result as a list of
[`arrow::Table`](https://arrow.apache.org/docs/r/reference/Table-class.html)s.
This is a low-level routine intended to be used by loaders for other
in-core formats, such as `Seurat`, which can be created from the
resulting Tables.

#### Usage

    SOMAExperimentAxisQuery$read(
      X_layers = NULL,
      obs_column_names = NULL,
      var_column_names = NULL
    )

#### Arguments

- `X_layers`:

  The name(s) of the `X` layer(s) to read and return.

- `obs_column_names, var_column_names`:

  Specify which column names in `var` and `obs` dataframes to read and
  return.

------------------------------------------------------------------------

### Method `to_sparse_matrix()`

Retrieve a collection layer as a sparse matrix with named dimensions.

Load any layer from the `X`, `obsm`, `varm`, `obsp`, or `varp`
collections as a [sparse
matrix](https://rdrr.io/pkg/Matrix/man/sparseMatrix-class.html).

By default the matrix dimensions are named using the `soma_joinid`
values in the specified layer's dimensions (e.g., `soma_dim_0`).
However, dimensions can be named using values from any `obs` or `var`
column that uniquely identifies each record by specifying the
`obs_index` and `var_index` arguments.

For layers in `obsm` or `varm`, the column axis (the axis not indexed by
“`obs`” or “`var`”) is set to the range of values present in
“`soma_dim_1`”; this ensures that gaps in this axis are preserved (eg.
when a query for “`obs`” that results in selecting entries that are all
zero for a given PC)

#### Usage

    SOMAExperimentAxisQuery$to_sparse_matrix(
      collection,
      layer_name,
      obs_index = NULL,
      var_index = NULL
    )

#### Arguments

- `collection`:

  The [`SOMACollection`](SOMACollection.md) containing the layer of
  interest, either: `"X"`, `"obsm"`, `"varm"`, `"obsp"`, or `"varp"`.

- `layer_name`:

  Name of the layer to retrieve from the `collection`.

- `obs_index, var_index`:

  Name of the column in `obs` or `var` (`var_index`) containing values
  that should be used as dimension labels in the resulting matrix.
  Whether the values are used as row or column labels depends on the
  selected `collection`:

  |            |                      |                      |
  |------------|----------------------|----------------------|
  | Collection | `obs_index`          | `var_index`          |
  | `X`        | row names            | column names         |
  | `obsm`     | row names            | ignored              |
  | `varm`     | ignored              | row names            |
  | `obsp`     | row and column names | ignored              |
  | `varp`     | ignored              | row and column names |

#### Returns

A
[`Matrix::sparseMatrix`](https://rdrr.io/pkg/Matrix/man/sparseMatrix-class.html)

------------------------------------------------------------------------

### Method `to_seurat()`

Loads the query as a
[`Seurat`](https://satijalab.github.io/seurat-object/reference/Seurat-class.html)
object

#### Usage

    SOMAExperimentAxisQuery$to_seurat(
      X_layers = c(counts = "counts", data = "logcounts"),
      obs_index = NULL,
      var_index = NULL,
      obs_column_names = NULL,
      var_column_names = NULL,
      obsm_layers = NULL,
      varm_layers = NULL,
      obsp_layers = NULL,
      drop_levels = FALSE,
      version = NULL
    )

#### Arguments

- `X_layers`:

  A named character of X layers to add to the Seurat assay where the
  names are the names of Seurat slots and the values are the names of
  layers within `X`; names should be one of:

  - “`counts`” to add the layer as `counts`

  - “`data`” to add the layer as `data`

  - “`scale.data`” to add the layer as `scale.data`

  At least one of “`counts`” or “`data`” is required

- `obs_index`:

  Name of column in `obs` to add as cell names; uses
  `paste0("cell", obs_joinids())` by default

- `var_index`:

  Name of column in `var` to add as feature names; uses
  `paste0("feature", var_joinids())` by default

- `obs_column_names`:

  Names of columns in `obs` to add as cell-level meta data; by default,
  loads all columns

- `var_column_names`:

  Names of columns in `var` to add as feature-level meta data; by
  default, loads all columns

- `obsm_layers`:

  Names of arrays in `obsm` to add as the cell embeddings; pass `FALSE`
  to suppress loading in any dimensional reductions; by default, loads
  all dimensional reduction information

- `varm_layers`:

  Named vector of arrays in `varm` to load in as the feature loadings;
  names must be names of arrays in `obsm` (eg.
  `varm_layers = c(X_pca = "PCs")`); pass `FALSE` to suppress loading in
  any feature loadings; will try to determine `varm_layers` from
  `obsm_layers`

- `obsp_layers`:

  Names of arrays in `obsp` to load in as
  [`Graph`](https://satijalab.github.io/seurat-object/reference/Graph-class.html)`s`;
  by default, loads all graphs

- `drop_levels`:

  Drop unused levels from `obs` and `var` factor columns

- `version`:

  Assay version to read query in as; by default, will try to infer assay
  type from the measurement itself

#### Returns

A
[`Seurat`](https://satijalab.github.io/seurat-object/reference/Seurat-class.html)
object

------------------------------------------------------------------------

### Method `to_seurat_assay()`

Loads the query as a Seurat
[`Assay`](https://satijalab.github.io/seurat-object/reference/Assay-class.html)

#### Usage

    SOMAExperimentAxisQuery$to_seurat_assay(
      X_layers = c(counts = "counts", data = "logcounts"),
      obs_index = NULL,
      var_index = NULL,
      var_column_names = NULL,
      drop_levels = FALSE,
      version = NULL
    )

#### Arguments

- `X_layers`:

  A named character of X layers to add to the Seurat assay where the
  names are the names of Seurat slots and the values are the names of
  layers within `X`; names should be one of:

  - “`counts`” to add the layer as `counts`

  - “`data`” to add the layer as `data`

  - “`scale.data`” to add the layer as `scale.data`

  At least one of “`counts`” or “`data`” is required

- `obs_index`:

  Name of column in `obs` to add as cell names; uses
  `paste0("cell", obs_joinids())` by default

- `var_index`:

  Name of column in `var` to add as feature names; uses
  `paste0("feature", var_joinids())` by default

- `var_column_names`:

  Names of columns in `var` to add as feature-level meta data; by
  default, loads all columns

- `drop_levels`:

  Drop unused levels from `var` factor columns

- `version`:

  Assay version to read query in as; by default, will try to infer assay
  type from the measurement itself

#### Returns

An
[`Assay`](https://satijalab.github.io/seurat-object/reference/Assay-class.html)
object

------------------------------------------------------------------------

### Method `to_seurat_reduction()`

Loads the query as a Seurat [dimensional
reduction](https://satijalab.github.io/seurat-object/reference/DimReduc-class.html)

#### Usage

    SOMAExperimentAxisQuery$to_seurat_reduction(
      obsm_layer,
      varm_layer = NULL,
      obs_index = NULL,
      var_index = NULL
    )

#### Arguments

- `obsm_layer`:

  Name of array in `obsm` to load as the cell embeddings

- `varm_layer`:

  Name of the array in `varm` to load as the feature loadings; by
  default, will try to determine `varm_layer` from `obsm_layer`

- `obs_index`:

  Name of column in `obs` to add as cell names; uses
  `paste0("cell", obs_joinids())` by default

- `var_index`:

  Name of column in `var` to add as feature names; uses
  `paste0("feature", var_joinids())` by default

#### Returns

A
[`DimReduc`](https://satijalab.github.io/seurat-object/reference/DimReduc-class.html)
object

------------------------------------------------------------------------

### Method `to_seurat_graph()`

Loads the query as a Seurat
[graph](https://satijalab.github.io/seurat-object/reference/Graph-class.html)

#### Usage

    SOMAExperimentAxisQuery$to_seurat_graph(obsp_layer, obs_index = NULL)

#### Arguments

- `obsp_layer`:

  Name of array in `obsp` to load as the graph

- `obs_index`:

  Name of column in `obs` to add as cell names; uses
  `paste0("cell", obs_joinids())` by default

#### Returns

A
[`Graph`](https://satijalab.github.io/seurat-object/reference/Graph-class.html)
object

------------------------------------------------------------------------

### Method `to_single_cell_experiment()`

Loads the query as a `SingleCellExperiment` object

#### Usage

    SOMAExperimentAxisQuery$to_single_cell_experiment(
      X_layers = NULL,
      obs_index = NULL,
      var_index = NULL,
      obs_column_names = NULL,
      var_column_names = NULL,
      obsm_layers = NULL,
      obsp_layers = NULL,
      varp_layers = NULL,
      drop_levels = FALSE
    )

#### Arguments

- `X_layers`:

  A character vector of X layers to add as assays in the main
  experiment; may optionally be named to set the name of the resulting
  assay (eg. `X_layers = c(counts = "raw")` will load in X layer “`raw`”
  as assay “`counts`”); by default, loads in all X layers

- `obs_index`:

  Name of column in `obs` to add as cell names; uses
  `paste0("cell", obs_joinids())` by default

- `var_index`:

  Name of column in `var` to add as feature names; uses
  `paste0("feature", var_joinids())` by default

- `obs_column_names`:

  Names of columns in `obs` to add as `colData`; by default, loads all
  columns

- `var_column_names`:

  Names of columns in `var` to add as `rowData`; by default, loads all
  columns

- `obsm_layers`:

  Names of arrays in `obsm` to add as the reduced dimensions; pass
  `FALSE` to suppress loading in any reduced dimensions; by default,
  loads all reduced dimensions

- `obsp_layers`:

  Names of arrays in `obsp` to load in as `SelfHits`; by default, loads
  all graphs

- `varp_layers`:

  Names of arrays in `varp` to load in as `SelfHits`; by default, loads
  all networks

- `drop_levels`:

  Drop unused levels from `obs` and `var` factor columns

#### Returns

A `SingleCellExperiment` object

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    SOMAExperimentAxisQuery$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
