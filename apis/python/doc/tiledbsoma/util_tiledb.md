<a id="tiledbsc.tiledbsoma/util_tiledb"></a>

# tiledbsc.tiledbsoma/util\_tiledb

<a id="tiledbsc.tiledbsoma/util_tiledb.tiledb_result_order_from_soma_result_order_non_indexed"></a>

#### tiledb\_result\_order\_from\_soma\_result\_order\_non\_indexed

```python
def tiledb_result_order_from_soma_result_order_non_indexed(
        soma_result_order: Optional[str]) -> Optional[str]
```

Maps SOMA-spec `result_order` syntax to TileDB-specific syntax, for non-indexed dataframes.

<a id="tiledbsc.tiledbsoma/util_tiledb.tiledb_result_order_from_soma_result_order_indexed"></a>

#### tiledb\_result\_order\_from\_soma\_result\_order\_indexed

```python
def tiledb_result_order_from_soma_result_order_indexed(
        soma_result_order: Optional[str]) -> Optional[str]
```

Maps SOMA-spec `result_order` syntax to TileDB-specific syntax, for indexed dataframes.

<a id="tiledbsc.tiledbsoma/util_tiledb.to_tiledb_supported_dtype"></a>

#### to\_tiledb\_supported\_dtype

```python
def to_tiledb_supported_dtype(dtype: np.dtype) -> np.dtype
```

A handful of types are cast into the TileDB type system.

<a id="tiledbsc.tiledbsoma/util_tiledb.to_tiledb_supported_array_type"></a>

#### to\_tiledb\_supported\_array\_type

```python
def to_tiledb_supported_array_type(x: T) -> T
```

Converts datatypes unrepresentable by TileDB into datatypes it can represent.
E.g., categorical strings -> string.

See also [https://docs.scipy.org/doc/numpy-1.10.1/reference/arrays.dtypes.html](https://docs.scipy.org/doc/numpy-1.10.1/reference/arrays.dtypes.html).

Preferentially converts to the underlying primitive type, as TileDB does not
support most complex types. NOTE: this does not support `datetime64` conversion.

Categoricals are a special case. If the underlying categorical type is a
primitive, convert to that. If the array contains NA/NaN (i.e. not in the
category, code == -1), raise error unless it is a float or string.

