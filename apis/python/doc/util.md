<a id="tiledbsc.util"></a>

# tiledbsc.util

<a id="tiledbsc.util.is_soma"></a>

#### is\_soma

```python
def is_soma(uri: str, ctx: Optional[tiledb.Ctx] = None) -> bool
```

Tells whether the URI points to a SOMA or not.

<a id="tiledbsc.util.is_soma_collection"></a>

#### is\_soma\_collection

```python
def is_soma_collection(uri: str, ctx: Optional[tiledb.Ctx] = None) -> bool
```

Tells whether the URI points to a SOMACollection or not.

<a id="tiledbsc.util.is_local_path"></a>

#### is\_local\_path

```python
def is_local_path(path: str) -> bool
```

Returns information about start time of an event. Nominally float seconds since the epoch,
but articulated here as being compatible with the format_elapsed function.

<a id="tiledbsc.util.get_start_stamp"></a>

#### get\_start\_stamp

```python
def get_start_stamp()
```

Returns information about start time of an event. Nominally float seconds since the epoch,
but articulated here as being compatible with the format_elapsed function.

<a id="tiledbsc.util.format_elapsed"></a>

#### format\_elapsed

```python
def format_elapsed(start_stamp, message: str)
```

Returns the message along with an elapsed-time indicator, with end time relative to start
start from `get_start_stamp`. Used for annotating elapsed time of a task.

<a id="tiledbsc.util.X_and_ids_to_sparse_matrix"></a>

#### X\_and\_ids\_to\_sparse\_matrix

```python
def X_and_ids_to_sparse_matrix(
    Xdf: pd.DataFrame,
    row_dim_name: str,
    col_dim_name: str,
    attr_name: str,
    row_labels: List[str],
    col_labels: List[str],
    return_as: str = "csr"
) -> Union[scipy.sparse.csr_matrix, scipy.sparse.csc_matrix]
```

This is needed when we read a TileDB X.df[:]. Since TileDB X is sparse 2D string-dimensioned,
the return value of which is a dict with three columns -- obs_id, var_id, and value. For
conversion to anndata, we need make a sparse COO/IJV-format array where the indices are
not strings but ints, matching the obs and var labels.
The `return_as` parameter must be one of `"csr"` or `"csc"`.

<a id="tiledbsc.util.triples_to_dense_df"></a>

#### triples\_to\_dense\_df

```python
def triples_to_dense_df(sparse_df: pd.DataFrame, fillna=0.0) -> pd.DataFrame
```

Output from X dataframe reads is in "triples" format, e.g. two index columns `obs_id` and `var_id`,
and data column `value`. This is the default format, and is appropriate for large, possibly sparse matrices.
However, sometimes we want a dense matrix with `obs_id` row labels, `var_id` column labels, and `value` data.
This function produces that.

<a id="tiledbsc.util.ETATracker"></a>

## ETATracker Objects

```python
class ETATracker()
```

Computes estimated time to completion for chunked writes.

<a id="tiledbsc.util.ETATracker.ingest_and_predict"></a>

#### ingest\_and\_predict

```python
def ingest_and_predict(chunk_percent: float, chunk_seconds: float) -> str
```

Updates from most recent chunk percent-done and chunk completion-seconds, then does a linear regression on all chunks done so far and estimates time to completion.

**Arguments**:

- `chunk_percent`: a percent done like 6.1 or 10.3.
- `chunk_seconds`: number of seconds it took to do the current chunk operation.

