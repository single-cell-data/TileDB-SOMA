<a id="tiledbsc.uns_array"></a>

# tiledbsc.uns\_array

<a id="tiledbsc.uns_array.UnsArray"></a>

## UnsArray Objects

```python
class UnsArray(TileDBArray)
```

Holds TileDB storage for an array obtained from the nested `anndata.uns` field.

<a id="tiledbsc.uns_array.UnsArray.__init__"></a>

#### \_\_init\_\_

```python
def __init__(uri: str, name: str, *, parent: Optional[TileDBGroup] = None)
```

See the TileDBObject constructor.

<a id="tiledbsc.uns_array.UnsArray.from_pandas_dataframe"></a>

#### from\_pandas\_dataframe

```python
def from_pandas_dataframe(df: pd.DataFrame) -> None
```

Ingests an `UnsArray` into TileDB storage, given a pandas.DataFrame.

<a id="tiledbsc.uns_array.UnsArray.from_numpy_ndarray"></a>

#### from\_numpy\_ndarray

```python
def from_numpy_ndarray(arr: np.ndarray) -> None
```

Writes a numpy.ndarray to a TileDB array, nominally for ingest of `uns` nested data from anndata
objects. Mostly tiledb.from_numpy, but with some necessary handling for data with UTF-8 values.

<a id="tiledbsc.uns_array.UnsArray.from_scipy_csr"></a>

#### from\_scipy\_csr

```python
def from_scipy_csr(csr: sp.csr_matrix) -> None
```

Convert ndarray/(csr|csc)matrix to coo_matrix and ingest into TileDB.

**Arguments**:

- `csr`: Matrix-like object coercible to a scipy coo_matrix.

<a id="tiledbsc.uns_array.UnsArray.create_empty_array_for_csr"></a>

#### create\_empty\_array\_for\_csr

```python
def create_empty_array_for_csr(attr_name: str, matrix_dtype: np.dtype,
                               nrows: int, ncols: int) -> None
```

Create a TileDB 2D sparse array with int dimensions and a single attribute.

Nominally used for uns data.

**Arguments**:

- `matrix_dtype`: datatype of the matrix
- `nrows`: number of rows in the matrix
- `ncols`: number of columns in the matrix

<a id="tiledbsc.uns_array.UnsArray.ingest_data_from_csr"></a>

#### ingest\_data\_from\_csr

```python
def ingest_data_from_csr(csr: sp.csr_matrix) -> None
```

Convert ndarray/(csr|csc)matrix to coo_matrix and ingest into TileDB.

**Arguments**:

- `csr`: Matrix-like object coercible to a scipy coo_matrix.

<a id="tiledbsc.uns_array.UnsArray.to_matrix"></a>

#### to\_matrix

```python
def to_matrix() -> np.ndarray
```

Reads an uns array from TileDB storage and returns a matrix -- currently, always as numpy.ndarray.

