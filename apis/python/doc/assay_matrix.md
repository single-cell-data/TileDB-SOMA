<a id="tiledbsc.assay_matrix"></a>

# tiledbsc.assay\_matrix

<a id="tiledbsc.assay_matrix.AssayMatrix"></a>

## AssayMatrix Objects

```python
class AssayMatrix(TileDBArray)
```

Wraps a TileDB sparse array with two string dimensions.
Used for `X`, `raw.X`, `obsp` elements, and `varp` elements.

<a id="tiledbsc.assay_matrix.AssayMatrix.__init__"></a>

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

See the TileDBObject constructor.

The `row_dataframe` and `col_dataframe` are nominally:

* `soma.obs` and `soma.var`, for `soma.X["data"]`
* `soma.obs` and `soma.raw.var`, for `soma.raw.X["data"]`
* `soma.obs` and `soma.obs`, for `soma.obsp` elements
* `soma.var` and `soma.var`, for `soma.obsp` elements

References to these objects are kept solely for obtaining dim labels for metadata
acquisition at runtime (e.g. shape). We retain references to these objects, rather
than taking in actualized ID-lists here in the constructor, for two reasons:

* We need to be able to set up a SOMA to write to, before it's been populated.
* For reading from an already-populated SOMA, we wish to avoid cache-coherency issues.

<a id="tiledbsc.assay_matrix.AssayMatrix.shape"></a>

#### shape

```python
def shape()
```

Returns a tuple with the number of rows and number of columns of the `AssayMatrix`.
In TileDB storage, these are string-indexed sparse arrays for which no `.shape()` exists,
but, we draw from the appropriate `obs`, `var`, `raw/var`, etc. as appropriate for a given matrix.

Note: currently implemented via data scan -- will be optimized for TileDB core 2.10.

<a id="tiledbsc.assay_matrix.AssayMatrix.dim_select"></a>

#### dim\_select

```python
def dim_select(obs_ids, var_ids)
```

Selects a slice out of the matrix with specified `obs_ids` and/or `var_ids`.
Either or both of the ID lists may be `None`, meaning, do not subselect along
that dimension. If both ID lists are `None`, the entire matrix is returned.

<a id="tiledbsc.assay_matrix.AssayMatrix.df"></a>

#### df

```python
def df(obs_ids=None, var_ids=None) -> pd.DataFrame
```

Keystroke-saving alias for `.dim_select()`. If either of `obs_ids` or `var_ids`
are provided, they're used to subselect; if not, the entire dataframe is returned.

<a id="tiledbsc.assay_matrix.AssayMatrix.from_matrix_and_dim_values"></a>

#### from\_matrix\_and\_dim\_values

```python
def from_matrix_and_dim_values(matrix, row_names, col_names) -> None
```

Imports a matrix -- nominally `scipy.sparse.csr_matrix` or `numpy.ndarray` -- into a TileDB
array which is used for `X`, `raw.X`, `obsp` members, and `varp` members.

<a id="tiledbsc.assay_matrix.AssayMatrix.ingest_data_whole"></a>

#### ingest\_data\_whole

```python
def ingest_data_whole(matrix, row_names, col_names) -> None
```

Convert `numpy.ndarray`, `scipy.sparse.csr_matrix`, or `scipy.sparse.csc_matrix` to COO matrix and ingest into TileDB.

**Arguments**:

- `matrix`: Matrix-like object coercible to a scipy COO matrix.
- `row_names`: List of row names.
- `col_names`: List of column names.

<a id="tiledbsc.assay_matrix.AssayMatrix.ingest_data_rows_chunked"></a>

#### ingest\_data\_rows\_chunked

```python
def ingest_data_rows_chunked(matrix, row_names, col_names) -> None
```

Convert csr_matrix to coo_matrix chunkwise and ingest into TileDB.

**Arguments**:

- `uri`: TileDB URI of the array to be written.
- `matrix`: csr_matrix.
- `row_names`: List of row names.
- `col_names`: List of column names.

<a id="tiledbsc.assay_matrix.AssayMatrix.ingest_data_cols_chunked"></a>

#### ingest\_data\_cols\_chunked

```python
def ingest_data_cols_chunked(matrix, row_names, col_names) -> None
```

Convert csc_matrix to coo_matrix chunkwise and ingest into TileDB.

**Arguments**:

- `uri`: TileDB URI of the array to be written.
- `matrix`: csc_matrix.
- `row_names`: List of row names.
- `col_names`: List of column names.

<a id="tiledbsc.assay_matrix.AssayMatrix.ingest_data_dense_rows_chunked"></a>

#### ingest\_data\_dense\_rows\_chunked

```python
def ingest_data_dense_rows_chunked(matrix, row_names, col_names) -> None
```

Convert dense matrix to coo_matrix chunkwise and ingest into TileDB.

**Arguments**:

- `uri`: TileDB URI of the array to be written.
- `matrix`: dense matrix.
- `row_names`: List of row names.
- `col_names`: List of column names.

<a id="tiledbsc.assay_matrix.AssayMatrix.to_csr_matrix"></a>

#### to\_csr\_matrix

```python
def to_csr_matrix(row_labels, col_labels)
```

Reads the TileDB array storage for the storage and returns a sparse CSR matrix.  The
row/columns labels should be `obs,var` labels if the `AssayMatrix` is `X`, or `obs,obs` labels if
the `AssayMatrix` is `obsp`, or `var,var` labels if the `AssayMatrix` is `varp`.
Note in all cases that TileDB will have sorted the row and column labels; they won't
be in the same order as they were in any anndata object which was used to create the
TileDB storage.

