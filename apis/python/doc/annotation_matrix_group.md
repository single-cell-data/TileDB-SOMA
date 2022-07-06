<a id="tiledbsc.annotation_matrix_group"></a>

# tiledbsc.annotation\_matrix\_group

<a id="tiledbsc.annotation_matrix_group.AnnotationMatrixGroup"></a>

## AnnotationMatrixGroup Objects

```python
class AnnotationMatrixGroup(TileDBGroup)
```

Nominally for soma obsm and varm. You can find element names using soma.obsm.keys(); you access
elements using soma.obsm['X_pca'] etc., or soma.obsm.X_pca if you prefer.  (The latter syntax is
possible when the element name doesn't have dashes, dots, etc. in it.)

<a id="tiledbsc.annotation_matrix_group.AnnotationMatrixGroup.__init__"></a>

#### \_\_init\_\_

```python
def __init__(uri: str, name: str, *, parent: Optional[TileDBGroup] = None)
```

See the TileDBObject constructor.

<a id="tiledbsc.annotation_matrix_group.AnnotationMatrixGroup.keys"></a>

#### keys

```python
def keys() -> List[str]
```

For `obsm` and `varm`, `.keys()` is a keystroke-saver for the more general group-member
accessor `._get_member_names()`.

<a id="tiledbsc.annotation_matrix_group.AnnotationMatrixGroup.__repr__"></a>

#### \_\_repr\_\_

```python
def __repr__() -> str
```

Default display of soma.obsm and soma.varm.

<a id="tiledbsc.annotation_matrix_group.AnnotationMatrixGroup.__iter__"></a>

#### \_\_iter\_\_

```python
def __iter__() -> Iterator[AnnotationMatrix]
```

Implements `for matrix in soma.obsm: ...` and `for matrix in soma.varm: ...`

<a id="tiledbsc.annotation_matrix_group.AnnotationMatrixGroup.__getattr__"></a>

#### \_\_getattr\_\_

```python
def __getattr__(name) -> Optional[AnnotationMatrix]
```

This is called on `soma.obsm.name` when `name` is not already an attribute.
This way you can do `soma.obsm.X_tsne` as an alias for `soma.obsm['X_tsne']`.

<a id="tiledbsc.annotation_matrix_group.AnnotationMatrixGroup.__getitem__"></a>

#### \_\_getitem\_\_

```python
def __getitem__(name) -> Optional[AnnotationMatrix]
```

Returns an `AnnotationMatrix` element at the given name within the group, or None if no such
member exists.  Overloads the `[...]` operator.

<a id="tiledbsc.annotation_matrix_group.AnnotationMatrixGroup.__contains__"></a>

#### \_\_contains\_\_

```python
def __contains__(name) -> bool
```

Implements the `in` operator, e.g. `"namegoeshere" in soma.obsm/soma.varm`.

<a id="tiledbsc.annotation_matrix_group.AnnotationMatrixGroup.add_matrix_from_matrix_and_dim_values"></a>

#### add\_matrix\_from\_matrix\_and\_dim\_values

```python
def add_matrix_from_matrix_and_dim_values(matrix, dim_values,
                                          matrix_name: str) -> None
```

Populates a component of the `obsm` or `varm` subgroup for a SOMA object.

**Arguments**:

- `matrix`: element of anndata.obsm, anndata.varm, or anndata.raw.varm.
- `dim_values`: anndata.obs_names, anndata.var_names, or anndata.raw.var_names.
- `matrix_name_name`: name of the matrix, like `"X_tsne"` or `"PCs"`.

<a id="tiledbsc.annotation_matrix_group.AnnotationMatrixGroup.remove"></a>

#### remove

```python
def remove(matrix_name: str) -> None
```

Removes a component of the `obsm` or `varm` subgroup for a SOMA object,
when invoked as `soma.obsm.remove("namegoeshere").

<a id="tiledbsc.annotation_matrix_group.AnnotationMatrixGroup.__delattr__"></a>

#### \_\_delattr\_\_

```python
def __delattr__(matrix_name: str) -> None
```

Removes a component of the `obsm` or `varm` subgroup for a SOMA object,
when invoked as `del soma.obsm.namegoeshere`.

<a id="tiledbsc.annotation_matrix_group.AnnotationMatrixGroup.__delitem__"></a>

#### \_\_delitem\_\_

```python
def __delitem__(matrix_name: str) -> None
```

Removes a component of the `obsm` or `varm` subgroup for a SOMA object,
when invoked as `del soma.obsm["namegoeshere"]`.

<a id="tiledbsc.annotation_matrix_group.AnnotationMatrixGroup.to_dict_of_csr"></a>

#### to\_dict\_of\_csr

```python
def to_dict_of_csr() -> Dict[str, scipy.sparse.csr_matrix]
```

Reads the obsm/varm group-member arrays into a dict from name to member array.
Member arrays are returned in sparse CSR format.

