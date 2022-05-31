<a id="tiledbsc.assay_matrix_group"></a>

# tiledbsc.assay\_matrix\_group

<a id="tiledbsc.assay_matrix_group.AssayMatrixGroup"></a>

## AssayMatrixGroup Objects

```python
class AssayMatrixGroup(TileDBGroup)
```

Nominally for `X`, `raw/X`, `obsp` elements, and `varp` elements.

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

<a id="tiledbsc.assay_matrix_group.AssayMatrixGroup.from_matrix_and_dim_values"></a>

#### from\_matrix\_and\_dim\_values

```python
def from_matrix_and_dim_values(matrix, row_names, col_names) -> None
```

Populates the `X` or `raw.X` subgroup for a `SOMA` object.  For `X` and `raw.X`, nominally `row_names` will be `anndata.obs_names` and `col_names` will be `anndata.var_names` or `anndata.raw.var_names`.  For `obsp` elements, both will be `anndata.obs_names`; for `varp elements, both will be `anndata.var_names`.

