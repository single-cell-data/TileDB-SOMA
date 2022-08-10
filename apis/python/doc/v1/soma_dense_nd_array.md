<a id="tiledbsc.v1/soma_dense_nd_array"></a>

# tiledbsc.v1/soma\_dense\_nd\_array

<a id="tiledbsc.v1/soma_dense_nd_array.SOMADenseNdArray"></a>

## SOMADenseNdArray Objects

```python
class SOMADenseNdArray(TileDBArray)
```

Represents ``X`` and others.

<a id="tiledbsc.v1/soma_dense_nd_array.SOMADenseNdArray.__init__"></a>

#### \_\_init\_\_

```python
def __init__(uri: str,
             *,
             name: Optional[str] = None,
             parent: Optional[SOMACollection] = None)
```

Also see the :class:`TileDBObject` constructor.

<a id="tiledbsc.v1/soma_dense_nd_array.SOMADenseNdArray.create"></a>

#### create

```python
def create(type: pa.DataType, shape: Union[NTuple, List[int]]) -> None
```

Create a SOMADenseNdArray named with the URI.

**Arguments**:

- `type`: an Arrow type defining the type of each element in the array. If the type is
unsupported, an error will be raised.
- `shape`: the length of each domain as a list, e.g., [100, 10]. All lengths must be in
the uint64 range.

<a id="tiledbsc.v1/soma_dense_nd_array.SOMADenseNdArray.__repr__"></a>

#### \_\_repr\_\_

```python
def __repr__() -> str
```

Default display of `SOMADenseNdArray`.

<a id="tiledbsc.v1/soma_dense_nd_array.SOMADenseNdArray.get_shape"></a>

#### get\_shape

```python
def get_shape() -> NTuple
```

Return length of each dimension, always a list of length ``ndims``

<a id="tiledbsc.v1/soma_dense_nd_array.SOMADenseNdArray.get_ndims"></a>

#### get\_ndims

```python
def get_ndims() -> int
```

Return number of index columns

<a id="tiledbsc.v1/soma_dense_nd_array.SOMADenseNdArray.get_is_sparse"></a>

#### get\_is\_sparse

```python
def get_is_sparse() -> bool
```

Returns ``False``.

<a id="tiledbsc.v1/soma_dense_nd_array.SOMADenseNdArray.read"></a>

#### read

```python
def read(slice: Any) -> Any
```

Read a user-specified subset of the object, and return as one or more Arrow.Tensor.

**Arguments**:

- `slice`: per-dimension slice, expressed as a scalar, a range, or a list of both.
- `partitions`: an optional [`SOMAReadPartitions`](`SOMAReadPartitions`) hint to indicate
how results should be organized.
- `result_order`: order of read results. Can be one of row-major or column-major.
The `read` operation will return a language-specific iterator over one or more Arrow Tensor
objects and information describing them, allowing the incremental processing of results larger
than available memory. The actual iterator used is delegated to language-specific SOMA specs. The
`DenseReadResult` should include:

* The coordinates of the slice (e.g., origin, shape)
* an Arrow.Tensor with the slice values

<a id="tiledbsc.v1/soma_dense_nd_array.SOMADenseNdArray.write"></a>

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

<a id="tiledbsc.v1/soma_dense_nd_array.SOMADenseNdArray.to_pandas"></a>

#### to\_pandas

```python
def to_pandas() -> pd.DataFrame
```

TODO: comment

<a id="tiledbsc.v1/soma_dense_nd_array.SOMADenseNdArray.from_matrix"></a>

#### from\_matrix

```python
def from_matrix(matrix: Matrix) -> None
```

Imports a matrix -- nominally `numpy.ndarray` -- into a TileDB
array which is used for `obsp` and `varp` matrices

