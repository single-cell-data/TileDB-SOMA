<a id="tiledbsc.tiledbsoma/soma_indexed_dataframe"></a>

# tiledbsc.tiledbsoma/soma\_indexed\_dataframe

<a id="tiledbsc.tiledbsoma/soma_indexed_dataframe.SOMAIndexedDataFrame"></a>

## SOMAIndexedDataFrame Objects

```python
class SOMAIndexedDataFrame(TileDBArray)
```

Represents ``obs``, ``var``, and others.

A SOMAIndexedDataFrame contains a "pseudo-column" called soma_rowid, of type uint64 and domain
[0,num_rows).  The soma_rowid pseudo-column contains a unique value for each row in the
SOMAIndexedDataFrame, and is intended to act as a join key for other objects, such as a SOMANdArray.

<a id="tiledbsc.tiledbsoma/soma_indexed_dataframe.SOMAIndexedDataFrame.__init__"></a>

#### \_\_init\_\_

```python
def __init__(uri: str,
             *,
             name: Optional[str] = None,
             parent: Optional[SOMACollection] = None)
```

See also the `TileDBOject` constructor.

<a id="tiledbsc.tiledbsoma/soma_indexed_dataframe.SOMAIndexedDataFrame.create"></a>

#### create

```python
def create(schema: pa.Schema,
           index_column_names: Optional[List[str]] = None) -> None
```

**Arguments**:

- `schema`: Arrow Schema defining the per-column schema. This schema must define all
columns, including columns to be named as index columns. The column name ``soma_rowid`` is
reserved for the pseudo-column of the same name. If the schema includes types unsupported by
the SOMA implementation, an error will be raised.
- `index_column_names`: A list of column names to use as user-defined index columns
(e.g., ``['cell_type', 'tissue_type']``). All named columns must exist in the schema, and at
least one index column name is required.

<a id="tiledbsc.tiledbsoma/soma_indexed_dataframe.SOMAIndexedDataFrame.__repr__"></a>

#### \_\_repr\_\_

```python
def __repr__() -> str
```

Default display of `SOMAIndexedDataFrame`.

<a id="tiledbsc.tiledbsoma/soma_indexed_dataframe.SOMAIndexedDataFrame.keys"></a>

#### keys

```python
def keys() -> List[str]
```

Returns the names of the columns when read back as a dataframe.
TODO: make it clear whether or not this will read back soma_rowid / soma_joinid.

<a id="tiledbsc.tiledbsoma/soma_indexed_dataframe.SOMAIndexedDataFrame.get_shape"></a>

#### get\_shape

```python
def get_shape() -> NTuple
```

Return length of each dimension, always a list of length ``ndims``

<a id="tiledbsc.tiledbsoma/soma_indexed_dataframe.SOMAIndexedDataFrame.get_ndims"></a>

#### get\_ndims

```python
def get_ndims() -> int
```

Return number of index columns

<a id="tiledbsc.tiledbsoma/soma_indexed_dataframe.SOMAIndexedDataFrame.get_index_column_names"></a>

#### get\_index\_column\_names

```python
def get_index_column_names() -> List[str]
```

Return index (dimension) column names.

<a id="tiledbsc.tiledbsoma/soma_indexed_dataframe.SOMAIndexedDataFrame.read"></a>

#### read

```python
def read(*,
         ids: Optional[Any] = None,
         value_filter: Optional[str] = None,
         column_names: Optional[Union[Sequence[str], str]] = None,
         result_order: Optional[str] = None) -> Iterator[pa.RecordBatch]
```

Read a user-defined subset of data, addressed by the dataframe indexing columns, optionally

filtered, and return results as one or more Arrow.RecordBatch.

**Arguments**:

- `ids`: for each index dimension, which rows to read. Defaults to 'all'.
- `column_names`: the named columns to read and return. Defaults to 'all'.
- `partitions`: an optional ``SOMAReadPartitions`` hint to indicate how results should be
organized.
- `result_order`: order of read results. This can be one of 'row-major', 'col-major', or
'unordered'.
- `value_filter`: an optional [value filter] to apply to the results. Defaults to no
filter.

**Indexing**: the `ids` parameter will support, per dimension: a list of values of the type
of the indexed column.

<a id="tiledbsc.tiledbsoma/soma_indexed_dataframe.SOMAIndexedDataFrame.read_all"></a>

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

<a id="tiledbsc.tiledbsoma/soma_indexed_dataframe.SOMAIndexedDataFrame.write"></a>

#### write

```python
def write(values: pa.RecordBatch) -> None
```

Write an Arrow.RecordBatch to the persistent object. As duplicate index values are not allowed,

index values already present in the object are overwritten and new index values are added.

**Arguments**:

- `values`: An Arrow.RecordBatch containing all columns, including the index columns. The
schema for the values must match the schema for the `SOMAIndexedDataFrame`.

<a id="tiledbsc.tiledbsoma/soma_indexed_dataframe.SOMAIndexedDataFrame.read_as_pandas"></a>

#### read\_as\_pandas

```python
def read_as_pandas(attrs: Optional[Sequence[str]] = None,
                   id_column_name: Optional[str] = None) -> Generator
```

For `to_anndata`, as well as for any interactive use where the user wants a Pandas dataframe.

<a id="tiledbsc.tiledbsoma/soma_indexed_dataframe.SOMAIndexedDataFrame.write_from_pandas"></a>

#### write\_from\_pandas

```python
def write_from_pandas(dataframe: pd.DataFrame,
                      index_column_names: List[str],
                      *,
                      extent: int = 2048,
                      id_column_name: Optional[str] = None) -> None
```

Populates the `obs` element of a SOMAExperiment object.

**Arguments**:

- `dataframe`: `anndata.obs`
- `extent`: TileDB `extent` parameter for the array schema.

