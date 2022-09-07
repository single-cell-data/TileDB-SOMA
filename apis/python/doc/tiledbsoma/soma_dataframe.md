<a id="tiledbsoma/soma_dataframe"></a>

# tiledbsoma/soma\_dataframe

<a id="tiledbsoma/soma_dataframe.SOMADataFrame"></a>

## SOMADataFrame Objects

```python
class SOMADataFrame(TileDBArray)
```

Represents ``obs``, ``var``, and others.

A `SOMADataFrame` contains a "pseudo-column" called `soma_rowid`, of type uint64 and domain
[0,num_rows).  The `soma_rowid` pseudo-column contains a unique value for each row in the
`SOMADataFrame`, and is intended to act as a join key for other objects, such as a `SOMASparseNdArray`.

<a id="tiledbsoma/soma_dataframe.SOMADataFrame.__init__"></a>

#### \_\_init\_\_

```python
def __init__(uri: str,
             *,
             name: Optional[str] = None,
             parent: Optional[SOMACollection] = None)
```

See also the `TileDBOject` constructor.

<a id="tiledbsoma/soma_dataframe.SOMADataFrame.create"></a>

#### create

```python
def create(schema: pa.Schema) -> None
```

**Arguments**:

- `schema`: Arrow Schema defining the per-column schema. This schema must define all
columns. The column name ``soma_rowid`` is reserved for the pseudo-column of the same name.
If the schema includes types unsupported by the SOMA implementation, an error will be
raised.

<a id="tiledbsoma/soma_dataframe.SOMADataFrame.__repr__"></a>

#### \_\_repr\_\_

```python
def __repr__() -> str
```

Default display of `SOMADataFrame`.

<a id="tiledbsoma/soma_dataframe.SOMADataFrame.keys"></a>

#### keys

```python
def keys() -> List[str]
```

Returns the names of the columns when read back as a dataframe.
TODO: make it clear whether or not this will read back `soma_rowid` / `soma_joinid`.

<a id="tiledbsoma/soma_dataframe.SOMADataFrame.get_shape"></a>

#### get\_shape

```python
def get_shape() -> NTuple
```

Return length of each dimension, always a list of length ``ndims``.

<a id="tiledbsoma/soma_dataframe.SOMADataFrame.get_ndims"></a>

#### get\_ndims

```python
def get_ndims() -> int
```

Return number of index columns.

<a id="tiledbsoma/soma_dataframe.SOMADataFrame.read"></a>

#### read

```python
def read(*,
         ids: Optional[Any] = None,
         value_filter: Optional[str] = None,
         column_names: Optional[Sequence[str]] = None,
         result_order: Optional[str] = None) -> Iterator[pa.RecordBatch]
```

Read a user-defined subset of data, addressed by the dataframe indexing column, optionally filtered, and return results as one or more `Arrow.RecordBatch`.

**Arguments**:

- `ids`: Which rows to read. Defaults to `None`, meaning no constraint -- all rows.
- `column_names`: the named columns to read and return. Defaults to `None`, meaning no constraint -- all column names.
- `partitions`: an optional ``SOMAReadPartitions`` hint to indicate how results should be
organized.
- `result_order`: order of read results.  This can be one of 'row-major', 'col-major', or
'unordered'.
- `value_filter`: an optional [value filter] to apply to the results. Defaults to no
filter.

**Indexing**: the `ids` parameter will support, per dimension: a row offset (uint), a
row-offset range (slice), or a list of both.

<a id="tiledbsoma/soma_dataframe.SOMADataFrame.read_all"></a>

#### read\_all

```python
def read_all(*,
             ids: Optional[Any] = None,
             value_filter: Optional[str] = None,
             column_names: Optional[Sequence[str]] = None,
             result_order: Optional[str] = None) -> pa.RecordBatch
```

This is a convenience method around `read`. It iterates the return value from `read`
and returns a concatenation of all the record batches found. Its nominal use is to
simply unit-test cases.

<a id="tiledbsoma/soma_dataframe.SOMADataFrame.write"></a>

#### write

```python
def write(values: pa.RecordBatch) -> None
```

Write an Arrow.RecordBatch to the persistent object.

**Arguments**:

- `values`: An Arrow.RecordBatch containing all columns, including the index columns. The
schema for the values must match the schema for the `SOMADataFrame`.

The ``values`` Arrow RecordBatch must contain a ``soma_rowid`` (uint64) column, indicating
which rows are being written.

<a id="tiledbsoma/soma_dataframe.SOMADataFrame.read_as_pandas"></a>

#### read\_as\_pandas

```python
def read_as_pandas(
        *,
        ids: Optional[Ids] = None,
        value_filter: Optional[str] = None,
        column_names: Optional[Sequence[str]] = None,
        result_order: Optional[str] = None,
        id_column_name: Optional[str] = None) -> Iterator[pd.DataFrame]
```

Reads from SOMA storage into memory.  For `to_anndata`, as well as for any interactive use
where the user wants a Pandas dataframe.  Returns a generator over dataframes for batched
read. See also `read_as_pandas_all` for a convenience wrapper.

TODO: params-list

<a id="tiledbsoma/soma_dataframe.SOMADataFrame.read_as_pandas_all"></a>

#### read\_as\_pandas\_all

```python
def read_as_pandas_all(*,
                       ids: Optional[Ids] = None,
                       value_filter: Optional[str] = None,
                       column_names: Optional[Sequence[str]] = None,
                       result_order: Optional[str] = None,
                       id_column_name: Optional[str] = None) -> pd.DataFrame
```

Reads from SOMA storage into memory.  Iterates over batches from `read_as_pandas`, concatenating
the output into a single dataframe.  Convenient for unit-test use; also, handy whenever
you're certain that the data being queried can be read entirely into memory.

<a id="tiledbsoma/soma_dataframe.SOMADataFrame.write_from_pandas"></a>

#### write\_from\_pandas

```python
def write_from_pandas(dataframe: pd.DataFrame,
                      *,
                      extent: int = 2048,
                      id_column_name: Optional[str] = None) -> None
```

Writes from memory to SOMA storage.

**Arguments**:

- `dataframe`: `anndata.obs` for example.
- `extent`: TileDB `extent` parameter for the array schema.

