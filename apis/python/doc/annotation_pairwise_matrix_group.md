<a id="tiledbsc.annotation_pairwise_matrix_group"></a>

# tiledbsc.annotation\_pairwise\_matrix\_group

<a id="tiledbsc.annotation_pairwise_matrix_group.AnnotationPairwiseMatrixGroup"></a>

## AnnotationPairwiseMatrixGroup Objects

```python
class AnnotationPairwiseMatrixGroup(TileDBGroup)
```

Nominally for soma obsp and varp. You can find element names using soma.obsp.keys(); you access
elements using soma.obsp['distances'] etc., or soma.obsp.distances if you prefer.  (The latter
syntax is possible when the element name doesn't have dashes, dots, etc. in it.)

<a id="tiledbsc.annotation_pairwise_matrix_group.AnnotationPairwiseMatrixGroup.__init__"></a>

#### \_\_init\_\_

```python
def __init__(uri: str,
             name: str,
             row_dataframe: AnnotationDataFrame,
             col_dataframe: AnnotationDataFrame,
             *,
             parent: Optional[TileDBGroup] = None)
```

See the `TileDBObject` constructor.
See `AssayMatrix` for the rationale behind retaining references to the `row_dataframe` and
`col_dataframe` objects.

<a id="tiledbsc.annotation_pairwise_matrix_group.AnnotationPairwiseMatrixGroup.keys"></a>

#### keys

```python
def keys() -> List[str]
```

For obsp and varp, `.keys()` is a keystroke-saver for the more general group-member
accessor `._get_member_names()`.

<a id="tiledbsc.annotation_pairwise_matrix_group.AnnotationPairwiseMatrixGroup.__repr__"></a>

#### \_\_repr\_\_

```python
def __repr__() -> str
```

Default display of soma.obsp and soma.varp.

<a id="tiledbsc.annotation_pairwise_matrix_group.AnnotationPairwiseMatrixGroup.__getattr__"></a>

#### \_\_getattr\_\_

```python
def __getattr__(name) -> Optional[AssayMatrix]
```

This is called on `soma.obsp.name` when `name` is not already an attribute.
This way you can do `soma.obsp.distances` as an alias for `soma.obsp['distances']`.

<a id="tiledbsc.annotation_pairwise_matrix_group.AnnotationPairwiseMatrixGroup.__getitem__"></a>

#### \_\_getitem\_\_

```python
def __getitem__(name) -> Optional[AssayMatrix]
```

Returns an `AssayMatrix` element at the given name within the group, or `None` if no such
member exists.  Overloads the `[...]` operator.

<a id="tiledbsc.annotation_pairwise_matrix_group.AnnotationPairwiseMatrixGroup.__contains__"></a>

#### \_\_contains\_\_

```python
def __contains__(name) -> bool
```

Implements `"namegoeshere" in soma.obsp/soma.varp`.

<a id="tiledbsc.annotation_pairwise_matrix_group.AnnotationPairwiseMatrixGroup.__iter__"></a>

#### \_\_iter\_\_

```python
def __iter__() -> Iterator[AssayMatrix]
```

Implements `for matrix in soma.obsp: ...` and `for matrix in soma.varp: ...`

<a id="tiledbsc.annotation_pairwise_matrix_group.AnnotationPairwiseMatrixGroup.remove"></a>

#### remove

```python
def remove(matrix_name: str) -> None
```

Removes a component of the `obsp` or `varp` subgroup for a SOMA object.
Implements `del soma.obsp['distances']` etc.

<a id="tiledbsc.annotation_pairwise_matrix_group.AnnotationPairwiseMatrixGroup.__delete__"></a>

#### \_\_delete\_\_

```python
def __delete__(matrix_name: str) -> None
```

Removes a component of the `obsp` or `varp` subgroup for a SOMA object.

<a id="tiledbsc.annotation_pairwise_matrix_group.AnnotationPairwiseMatrixGroup.add_matrix_from_matrix_and_dim_values"></a>

#### add\_matrix\_from\_matrix\_and\_dim\_values

```python
def add_matrix_from_matrix_and_dim_values(matrix, dim_values,
                                          matrix_name: str) -> None
```

Populates a component of the `obsp` or `varp` subgroup for a SOMA object.

**Arguments**:

- `matrix`: element of anndata.obsp or anndata.varp.
- `dim_values`: anndata.obs_names or anndata.var_names.
- `matrix_name_name`: name of the matrix, like `"distances"`.

<a id="tiledbsc.annotation_pairwise_matrix_group.AnnotationPairwiseMatrixGroup.to_dict_of_csr"></a>

#### to\_dict\_of\_csr

```python
def to_dict_of_csr(obs_df_index,
                   var_df_index) -> Dict[str, scipy.sparse.csr_matrix]
```

Reads the `obsp` or `varp` group-member arrays into a dict from name to member array.
Member arrays are returned in sparse CSR format.

