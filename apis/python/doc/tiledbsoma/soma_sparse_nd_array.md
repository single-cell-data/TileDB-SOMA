<a id="tiledbsoma/soma_sparse_nd_array"></a>

# tiledbsoma/soma\_sparse\_nd\_array

<a id="tiledbsoma/soma_sparse_nd_array.SOMASparseNdArray"></a>

## SOMASparseNdArray Objects

```python
class SOMASparseNdArray(TileDBArray)
```

Represents ``X`` and others.

<a id="tiledbsoma/soma_sparse_nd_array.SOMASparseNdArray.__init__"></a>

#### \_\_init\_\_

```python
def __init__(uri: str,
             *,
             name: Optional[str] = None,
             parent: Optional[SOMACollection] = None)
```

Also see the `TileDBObject` constructor.

<a id="tiledbsoma/soma_sparse_nd_array.SOMASparseNdArray.create"></a>

#### create

```python
def create(type: pa.DataType, shape: Union[NTuple, List[int]]) -> None
```

Create a `SOMASparseNdArray` named with the URI.

**Arguments**:

- `type`: an Arrow type defining the type of each element in the array. If the type is
unsupported, an error will be raised.
- `shape`: the length of each domain as a list, e.g., [100, 10]. All lengths must be in
the uint64 range.

<a id="tiledbsoma/soma_sparse_nd_array.SOMASparseNdArray.__repr__"></a>

#### \_\_repr\_\_

```python
def __repr__() -> str
```

Default display of `SOMASparseNdArray`.

<a id="tiledbsoma/soma_sparse_nd_array.SOMASparseNdArray.get_shape"></a>

#### get\_shape

```python
def get_shape() -> NTuple
```

Return length of each dimension, always a list of length ``ndims``

<a id="tiledbsoma/soma_sparse_nd_array.SOMASparseNdArray.get_ndims"></a>

#### get\_ndims

```python
def get_ndims() -> int
```

Return number of index columns

<a id="tiledbsoma/soma_sparse_nd_array.SOMASparseNdArray.get_is_sparse"></a>

#### get\_is\_sparse

```python
def get_is_sparse() -> bool
```

Returns ``True``.

<a id="tiledbsoma/soma_sparse_nd_array.SOMASparseNdArray.read"></a>

#### read

```python
def read(*,
         row_ids: Optional[Sequence[int]] = None,
         col_ids: Optional[Sequence[int]] = None,
         result_order: Optional[str] = None) -> Iterator[pa.RecordBatch]
```

TODO: comment

<a id="tiledbsoma/soma_sparse_nd_array.SOMASparseNdArray.read_as_pandas"></a>

#### read\_as\_pandas

```python
def read_as_pandas(*,
                   row_ids: Optional[Sequence[int]] = None,
                   col_ids: Optional[Sequence[int]] = None,
                   set_index: Optional[bool] = False) -> pd.DataFrame
```

TODO: comment

<a id="tiledbsoma/soma_sparse_nd_array.SOMASparseNdArray.read_all"></a>

#### read\_all

```python
def read_all(*,
             row_ids: Optional[Sequence[int]] = None,
             col_ids: Optional[Sequence[int]] = None,
             result_order: Optional[str] = None) -> pa.RecordBatch
```

This is a convenience method around `read`. It iterates the return value from `read`
and returns a concatenation of all the record batches found. Its nominal use is to
simply unit-test cases.

<a id="tiledbsoma/soma_sparse_nd_array.SOMASparseNdArray.read_as_pandas_all"></a>

#### read\_as\_pandas\_all

```python
def read_as_pandas_all(*,
                       row_ids: Optional[Sequence[int]] = None,
                       col_ids: Optional[Sequence[int]] = None,
                       set_index: Optional[bool] = False) -> pa.RecordBatch
```

This is a convenience method around `read_as_pandas`. It iterates the return value from
`read_as_pandas` and returns a concatenation of all the record batches found. Its nominal
use is to simply unit-test cases.

<a id="tiledbsoma/soma_sparse_nd_array.SOMASparseNdArray.write"></a>

#### write

```python
def write(tensor: Union[pa.SparseCOOTensor, pa.SparseCSFTensor]) -> None
```

Write an `Arrow.Tensor` to the persistent object. As duplicate index values are not allowed, index

values already present in the object are overwritten and new index values are added.

**Arguments**:

- `values`: an `Arrow.SparseTensor` containing values to be written. The type of elements in `values`
must match the type of the `SOMASparseNdArray`.

<a id="tiledbsoma/soma_sparse_nd_array.SOMASparseNdArray.from_matrix"></a>

#### from\_matrix

```python
def from_matrix(matrix: Matrix) -> None
```

Imports a matrix -- nominally `scipy.sparse.csr_matrix` or `numpy.ndarray` -- into a TileDB
array which is used for `X`, `obsm`, and `varm` matrices

