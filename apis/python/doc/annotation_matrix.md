<a id="tiledbsc.annotation_matrix"></a>

# tiledbsc.annotation\_matrix

<a id="tiledbsc.annotation_matrix.AnnotationMatrix"></a>

## AnnotationMatrix Objects

```python
class AnnotationMatrix(TileDBArray)
```

Nominally for obsm and varm group elements within a soma.

<a id="tiledbsc.annotation_matrix.AnnotationMatrix.__init__"></a>

#### \_\_init\_\_

```python
def __init__(uri: str,
             name: str,
             dim_name: str,
             *,
             parent: Optional[TileDBGroup] = None)
```

See the TileDBObject constructor.

<a id="tiledbsc.annotation_matrix.AnnotationMatrix.shape"></a>

#### shape

```python
def shape() -> Tuple[int, int]
```

Returns a tuple with the number of rows and number of columns of the `AnnotationMatrix`.
The row-count is the number of obs_ids (for `obsm` elements) or the number of var_ids (for
`varm` elements).  The column-count is the number of columns/attributes in the dataframe.

Note: currently implemented via data scan -- will be optimized for TileDB core 2.10.

<a id="tiledbsc.annotation_matrix.AnnotationMatrix.dim_select"></a>

#### dim\_select

```python
def dim_select(ids: Optional[Ids]) -> pd.DataFrame
```

Selects a slice out of the array with specified `obs_ids` (for `obsm` elements) or
`var_ids` (for `varm` elements).  If `ids` is `None`, the entire array is returned.

<a id="tiledbsc.annotation_matrix.AnnotationMatrix.df"></a>

#### df

```python
def df(ids: Optional[Ids] = None) -> pd.DataFrame
```

Keystroke-saving alias for `.dim_select()`. If `ids` are provided, they're used
to subselect; if not, the entire dataframe is returned.

<a id="tiledbsc.annotation_matrix.AnnotationMatrix.from_matrix_and_dim_values"></a>

#### from\_matrix\_and\_dim\_values

```python
def from_matrix_and_dim_values(matrix: Union[pd.DataFrame, Matrix],
                               dim_values: Labels) -> None
```

Populates an array in the obsm/ or varm/ subgroup for a SOMA object.

**Arguments**:

- `matrix`: anndata.obsm['foo'], anndata.varm['foo'], or anndata.raw.varm['foo'].
- `dim_values`: anndata.obs_names, anndata.var_names, or anndata.raw.var_names.

