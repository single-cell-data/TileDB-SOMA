<a id="tiledbsc.annotation_dataframe"></a>

# tiledbsc.annotation\_dataframe

<a id="tiledbsc.annotation_dataframe.AnnotationDataFrame"></a>

## AnnotationDataFrame Objects

```python
class AnnotationDataFrame(TileDBArray)
```

Nominally for `obs` and `var` data within a soma. These have one string dimension, and multiple attributes.

<a id="tiledbsc.annotation_dataframe.AnnotationDataFrame.__init__"></a>

#### \_\_init\_\_

```python
def __init__(uri: str, name: str, parent: Optional[TileDBGroup] = None)
```

See the TileDBObject constructor.

<a id="tiledbsc.annotation_dataframe.AnnotationDataFrame.shape"></a>

#### shape

```python
def shape()
```

Returns a tuple with the number of rows and number of columns of the `AnnotationDataFrame`.
The row-count is the number of obs_ids (for `obs`) or the number of var_ids (for `var`).
The column-count is the number of columns/attributes in the dataframe.

<a id="tiledbsc.annotation_dataframe.AnnotationDataFrame.ids"></a>

#### ids

```python
def ids() -> List[str]
```

Returns the `obs_ids` in the matrix (for `obs`) or the `var_ids` (for `var`).

<a id="tiledbsc.annotation_dataframe.AnnotationDataFrame.keys"></a>

#### keys

```python
def keys() -> List[str]
```

Returns the column names for the `obs` or `var` dataframe.  For obs and varp, `.keys()` is a
keystroke-saver for the more general array-schema accessor `attr_names`.

<a id="tiledbsc.annotation_dataframe.AnnotationDataFrame.keyset"></a>

#### keyset

```python
def keyset() -> Set[str]
```

Same as `.keys` but returns as set.

<a id="tiledbsc.annotation_dataframe.AnnotationDataFrame.dim_select"></a>

#### dim\_select

```python
def dim_select(ids, attrs=None)
```

Selects a slice out of the dataframe with specified `obs_ids` (for `obs`) or `var_ids` (for
`var`).  If `ids` is `None`, the entire dataframe is returned.  Similarly, if `attrs` are
provided, they're used for the query; else, all attributes are returned.

<a id="tiledbsc.annotation_dataframe.AnnotationDataFrame.df"></a>

#### df

```python
def df(ids=None, attrs=None) -> pd.DataFrame
```

Keystroke-saving alias for `.dim_select()`. If `ids` are provided, they're used
to subselect; if not, the entire dataframe is returned. If `attrs` are provided,
they're used for the query; else, all attributes are returned.

<a id="tiledbsc.annotation_dataframe.AnnotationDataFrame.query"></a>

#### query

```python
def query(query_string, ids=None, attrs=None)
```

Selects from obs/var using a TileDB-Py `QueryCondition` string such as `cell_type ==
"blood"`.  If `attrs` is `None`, returns all column names in the dataframe; use `[]` for
`attrs` to select none of them.  Any column names specified in the `query_string` must be
included in `attrs` if `attrs` is not `None`.  Returns `None` if the slice is empty.

<a id="tiledbsc.annotation_dataframe.AnnotationDataFrame.from_dataframe"></a>

#### from\_dataframe

```python
def from_dataframe(dataframe: pd.DataFrame, extent: int = 2048) -> None
```

Populates the `obs` or `var` subgroup for a SOMA object.

**Arguments**:

- `dataframe`: `anndata.obs`, `anndata.var`, `anndata.raw.var`.
- `extent`: TileDB `extent` parameter for the array schema.

