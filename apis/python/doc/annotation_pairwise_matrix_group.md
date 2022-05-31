<a id="tiledbsc.annotation_pairwise_matrix_group"></a>

# tiledbsc.annotation\_pairwise\_matrix\_group

<a id="tiledbsc.annotation_pairwise_matrix_group.AnnotationPairwiseMatrixGroup"></a>

## AnnotationPairwiseMatrixGroup Objects

```python
class AnnotationPairwiseMatrixGroup(TileDBGroup)
```

Nominally for soma obsp and varp.

<a id="tiledbsc.annotation_pairwise_matrix_group.AnnotationPairwiseMatrixGroup.__init__"></a>

#### \_\_init\_\_

```python
def __init__(uri: str,
             name: str,
             row_dataframe: AnnotationDataFrame,
             col_dataframe: AnnotationDataFrame,
             parent: Optional[TileDBGroup] = None)
```

See the `TileDBObject` constructor.
See `AssayMatrix` for the rationale behind retaining references to the `row_dataframe` and
`col_dataframe` objects.

<a id="tiledbsc.annotation_pairwise_matrix_group.AnnotationPairwiseMatrixGroup.keys"></a>

#### keys

```python
def keys()
```

For obsp and varp, `.keys()` is a keystroke-saver for the more general group-member
accessor `._get_member_names()`.

<a id="tiledbsc.annotation_pairwise_matrix_group.AnnotationPairwiseMatrixGroup.__iter__"></a>

#### \_\_iter\_\_

```python
def __iter__() -> List[AssayMatrix]
```

Implements `for matrix in soma.obsp: ...` and `for matrix in soma.varp: ...`

<a id="tiledbsc.annotation_pairwise_matrix_group.AnnotationPairwiseMatrixGroup.from_matrices_and_dim_values"></a>

#### from\_matrices\_and\_dim\_values

```python
def from_matrices_and_dim_values(annotation_pairwise_matrices, dim_values)
```

Populates the `obsp` or `varp` subgroup for a SOMA object, then writes all the components

arrays under that group.

**Arguments**:

- `annotation_pairwise_matrices`: anndata.obsp, anndata.varp, or anndata.raw.varp.
- `dim_values`: anndata.obs_names, anndata.var_names, or anndata.raw.var_names.

<a id="tiledbsc.annotation_pairwise_matrix_group.AnnotationPairwiseMatrixGroup.to_dict_of_csr"></a>

#### to\_dict\_of\_csr

```python
def to_dict_of_csr() -> Dict[str, scipy.sparse.csr_matrix]
```

Reads the `obsm` or `varm` group-member arrays into a dict from name to member array.
Member arrays are returned in sparse CSR format.

<a id="tiledbsc.annotation_pairwise_matrix_group.AnnotationPairwiseMatrixGroup.__getitem__"></a>

#### \_\_getitem\_\_

```python
def __getitem__(name)
```

Returns an `AssayMatrix` element at the given name within the group, or `None` if no such
member exists.  Overloads the `[...]` operator.

<a id="tiledbsc.annotation_pairwise_matrix_group.AnnotationPairwiseMatrixGroup.__contains__"></a>

#### \_\_contains\_\_

```python
def __contains__(name)
```

Implements `"namegoeshere" in soma.obsp/soma.varp`.

