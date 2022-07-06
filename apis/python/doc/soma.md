<a id="tiledbsc.soma"></a>

# tiledbsc.soma

<a id="tiledbsc.soma.SOMA"></a>

## SOMA Objects

```python
class SOMA(TileDBGroup)
```

Single-cell group
Class for representing a group of TileDB groups/arrays that constitute an SOMA ('stack of matrices, annotated')
which includes:

* `X` (group of `AssayMatrixGroup`): a group of one or more labeled 2D sparse arrays that share the same dimensions.
* `obs` (`AnnotationDataframe`): 1D labeled array with column labels for `X`
* `var` (`AnnotationDataframe`): 1D labeled array with row labels for `X`
* `obsm` (group of `AnnotationMatrix`): multi-attribute arrays keyed by IDs of `obs`
* `varm` (group of `AnnotationMatrix`): multi-attribute arrays keyed by IDs of `var`
* `obsp` (group of `AnnotationMatrix`): 2D arrays keyed by IDs of `obs`
* `varp` (group of `AnnotationMatrix`): 2D arrays keyed by IDs of `var`
* `raw`: contains raw versions of `X` and `varm`
* `uns`: nested, unstructured data

Convenience accessors include:

* `soma.obs_keys()` for `soma.obs_names` for `soma.obs.ids()`
* `soma.var_keys()` for `soma.var_names` for `soma.var.ids()`
* `soma.n_obs` for `soma.obs.shape()[0]`
* `soma.n_var` for `soma.var.shape()[0]`

<a id="tiledbsc.soma.SOMA.__init__"></a>

#### \_\_init\_\_

```python
def __init__(uri: str,
             *,
             name=None,
             soma_options: Optional[SOMAOptions] = None,
             config: Optional[tiledb.Config] = None,
             ctx: Optional[tiledb.Ctx] = None,
             parent: Optional[TileDBGroup] = None)
```

Create a new SOMA object. The existing array group is opened at the specified array `uri` if one is present, otherwise a new array group is created.

**Arguments**:

- `uri`: URI of the TileDB group

<a id="tiledbsc.soma.SOMA.__repr__"></a>

#### \_\_repr\_\_

```python
def __repr__() -> str
```

Default display of SOMA.

<a id="tiledbsc.soma.SOMA.__getattr__"></a>

#### \_\_getattr\_\_

```python
def __getattr__(name)
```

This is called on `soma.name` when `name` is not already an attribute.
This is used for `soma.n_obs`, etc.

<a id="tiledbsc.soma.SOMA.obs_keys"></a>

#### obs\_keys

```python
def obs_keys() -> List[str]
```

An alias for `soma.obs.ids()`.

<a id="tiledbsc.soma.SOMA.var_keys"></a>

#### var\_keys

```python
def var_keys() -> List[str]
```

An alias for `soma.var.ids()`.

<a id="tiledbsc.soma.SOMA.get_obs_value_counts"></a>

#### get\_obs\_value\_counts

```python
def get_obs_value_counts(obs_label: str) -> pd.DataFrame
```

Given an obs label, e.g. `cell_type`, returns a dataframe count the number of different
values for that label in the SOMA.

<a id="tiledbsc.soma.SOMA.get_var_value_counts"></a>

#### get\_var\_value\_counts

```python
def get_var_value_counts(var_label: str) -> pd.DataFrame
```

Given an var label, e.g. `feature_name`, returns a dataframe count the number of different
values for that label in the SOMA.

<a id="tiledbsc.soma.SOMA.dim_slice"></a>

#### dim\_slice

```python
def dim_slice(obs_ids, var_ids) -> Optional[SOMASlice]
```

Subselects the SOMA's obs, var, and X/data using the specified obs_ids and var_ids.
Using a value of `None` for obs_ids means use all obs_ids, and likewise for var_ids.
Returns `None` for empty slice.

<a id="tiledbsc.soma.SOMA.query"></a>

#### query

```python
def query(*,
          obs_attrs: Optional[List[str]] = None,
          obs_query_string: Optional[str] = None,
          obs_ids: Optional[List[str]] = None,
          var_attrs: Optional[List[str]] = None,
          var_query_string: Optional[str] = None,
          var_ids: Optional[List[str]] = None) -> Optional[SOMASlice]
```

Subselects the SOMA's obs, var, and X/data using the specified queries on obs and var.
Queries use the TileDB-Py `QueryCondition` API.

If `obs_query_string` is `None`, the `obs` dimension is not filtered and all of `obs` is
used; similiarly for `var`.

If `obs_attrs` or `var_attrs` are unspecified, the slice will take all `obs`/`var` attributes
from the source SOMAs; if they are specified, the slice will take the specified `obs`/`var`

<a id="tiledbsc.soma.SOMA.from_soma_slice"></a>

#### from\_soma\_slice

```python
@classmethod
def from_soma_slice(cls,
                    soma_slice: SOMASlice,
                    uri: str,
                    name=None,
                    soma_options: Optional[SOMAOptions] = None,
                    config: Optional[tiledb.Config] = None,
                    ctx: Optional[tiledb.Ctx] = None,
                    parent: Optional[TileDBGroup] = None) -> SOMA
```

Constructs `SOMA` storage from a given in-memory `SOMASlice` object.

