<a id="tiledbsoma/soma_dense_nd_array"></a>

# tiledbsoma/soma\_dense\_nd\_array

<a id="tiledbsoma/soma_dense_nd_array.SOMADenseNdArray"></a>

## SOMADenseNdArray Objects

```python
class SOMADenseNdArray(TileDBArray)
```

Represents ``X`` and others.

<a id="tiledbsoma/soma_dense_nd_array.SOMADenseNdArray.__init__"></a>

#### \_\_init\_\_

```python
def __init__(uri: str,
             *,
             name: Optional[str] = None,
             parent: Optional[SOMACollection] = None)
```

Also see the `TileDBObject` constructor.

<a id="tiledbsoma/soma_dense_nd_array.SOMADenseNdArray.create"></a>

#### create

```python
def create(type: pa.DataType, shape: Union[NTuple, List[int]]) -> None
```

Create a `SOMADenseNdArray` named with the URI.

**Arguments**:

- `type`: an Arrow type defining the type of each element in the array. If the type is
unsupported, an error will be raised.
- `shape`: the length of each domain as a list, e.g., [100, 10]. All lengths must be in
the uint64 range.

<a id="tiledbsoma/soma_dense_nd_array.SOMADenseNdArray.__repr__"></a>

#### \_\_repr\_\_

```python
def __repr__() -> str
```

Default display of `SOMADenseNdArray`.

<a id="tiledbsoma/soma_dense_nd_array.SOMADenseNdArray.get_shape"></a>

#### get\_shape

```python
def get_shape() -> NTuple
```

Return length of each dimension, always a list of length ``ndims``

<a id="tiledbsoma/soma_dense_nd_array.SOMADenseNdArray.get_ndims"></a>

#### get\_ndims

```python
def get_ndims() -> int
```

Return number of index columns

<a id="tiledbsoma/soma_dense_nd_array.SOMADenseNdArray.get_is_sparse"></a>

#### get\_is\_sparse

```python
def get_is_sparse() -> bool
```

Returns ``False``.

<a id="tiledbsoma/soma_dense_nd_array.SOMADenseNdArray.read"></a>

#### read

```python
def read(row_ids: Optional[Sequence[int]] = None,
         col_ids: Optional[Sequence[int]] = None,
         result_order: Optional[str] = None) -> Any
```

Read a user-specified subset of the object, and return as one or more Arrow.Tensor.

**Arguments**:

- `ids`: per-dimension slice, expressed as a scalar, a range, or a list of both.
- `partitions`: an optional [`SOMAReadPartitions`](`SOMAReadPartitions`) hint to indicate
how results should be organized.
- `result_order`: order of read results. Can be one of row-major or column-major.
The `read` operation will return a language-specific iterator over one or more Arrow Tensor
objects and information describing them, allowing the incremental processing of results larger
than available memory. The actual iterator used is delegated to language-specific SOMA specs. The
`DenseReadResult` should include:

* The coordinates of the slice (e.g., origin, shape)
* an Arrow.Tensor with the slice values

<a id="tiledbsoma/soma_dense_nd_array.SOMADenseNdArray.read_as_pandas"></a>

#### read\_as\_pandas

```python
def read_as_pandas(*,
                   row_ids: Optional[Sequence[int]] = None,
                   col_ids: Optional[Sequence[int]] = None,
                   set_index: Optional[bool] = False) -> pd.DataFrame
```

TODO: comment

<a id="tiledbsoma/soma_dense_nd_array.SOMADenseNdArray.read_all"></a>

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

<a id="tiledbsoma/soma_dense_nd_array.SOMADenseNdArray.read_as_pandas_all"></a>

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

<a id="tiledbsoma/soma_dense_nd_array.SOMADenseNdArray.write"></a>

#### write

```python
def write(coords: Any, values: pa.Tensor) -> None
```

Write an Arrow.Tensor to the persistent object. As duplicate index values are not allowed, index

values already present in the object are overwritten and new index values are added.

**Arguments**:

- `coords`: location at which to write the tensor
- `values`: an Arrow.Tensor containing values to be written. The type of elements in `values` must
match the type of the SOMADenseNdArray.

<a id="tiledbsoma/soma_dense_nd_array.SOMADenseNdArray.from_matrix"></a>

#### from\_matrix

```python
def from_matrix(matrix: Matrix) -> None
```

Imports a matrix -- nominally `numpy.ndarray` -- into a TileDB
array which is used for `obsp` and `varp` matrices

