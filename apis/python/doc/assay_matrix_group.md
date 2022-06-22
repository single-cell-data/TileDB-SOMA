<a id="tiledbsc.assay_matrix_group"></a>

# tiledbsc.assay\_matrix\_group

<a id="tiledbsc.assay_matrix_group.AssayMatrixGroup"></a>

## AssayMatrixGroup Objects

```python
class AssayMatrixGroup(TileDBGroup)
```

Nominally for `X` and `raw/X` elements.  You can find element names using soma.X.keys(); you
access elements using soma.X['data'] etc., or soma.X.data if you prefer.  (The latter syntax is
possible when the element name doesn't have dashes, dots, etc. in it.)

<a id="tiledbsc.assay_matrix_group.AssayMatrixGroup.__init__"></a>

#### \_\_init\_\_

```python
def __init__(uri: str,
             name: str,
             row_dim_name: str,
             col_dim_name: str,
             row_dataframe: AnnotationDataFrame,
             col_dataframe: AnnotationDataFrame,
             parent: Optional[TileDBGroup] = None)
```

See the `TileDBObject` constructor.

See `AssayMatrix` for the rationale behind retaining references to the `row_dataframe` and
`col_dataframe` objects.

<a id="tiledbsc.assay_matrix_group.AssayMatrixGroup.keys"></a>

#### keys

```python
def keys() -> List[str]
```

For `obsm` and `varm`, `.keys()` is a keystroke-saver for the more general group-member
accessor `._get_member_names()`.

<a id="tiledbsc.assay_matrix_group.AssayMatrixGroup.__getattr__"></a>

#### \_\_getattr\_\_

```python
def __getattr__(name) -> AssayMatrix
```

This is called on `soma.X.name` when `name` is not already an attribute.
This way you can do `soma.X.data` as an alias for `soma.X['data']`.

<a id="tiledbsc.assay_matrix_group.AssayMatrixGroup.__iter__"></a>

#### \_\_iter\_\_

```python
def __iter__() -> List[AssayMatrix]
```

Implements `for matrix in soma.obsm: ...` and `for matrix in soma.varm: ...`

<a id="tiledbsc.assay_matrix_group.AssayMatrixGroup.__getitem__"></a>

#### \_\_getitem\_\_

```python
def __getitem__(name) -> AssayMatrix
```

Returns an `AnnotationMatrix` element at the given name within the group, or None if no such
member exists.  Overloads the `[...]` operator.

<a id="tiledbsc.assay_matrix_group.AssayMatrixGroup.__contains__"></a>

#### \_\_contains\_\_

```python
def __contains__(name) -> bool
```

Implements the `in` operator, e.g. `"data" in soma.X`.

<a id="tiledbsc.assay_matrix_group.AssayMatrixGroup.add_layer_from_matrix_and_dim_values"></a>

#### add\_layer\_from\_matrix\_and\_dim\_values

```python
def add_layer_from_matrix_and_dim_values(matrix,
                                         row_names: str,
                                         col_names: str,
                                         layer_name="data") -> None
```

Populates the `X` or `raw.X` subgroup for a `SOMA` object.  For `X` and `raw.X`, nominally `row_names` will be `anndata.obs_names` and `col_names` will be `anndata.var_names` or `anndata.raw.var_names`.  For `obsp` elements, both will be `anndata.obs_names`; for `varp elements, both will be `anndata.var_names`.

