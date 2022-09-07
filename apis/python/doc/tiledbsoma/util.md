<a id="tiledbsoma/util"></a>

# tiledbsoma/util

<a id="tiledbsoma/util.get_start_stamp"></a>

#### get\_start\_stamp

```python
def get_start_stamp() -> float
```

Returns information about start time of an event. Nominally float seconds since the epoch,
but articulated here as being compatible with the format_elapsed function.

<a id="tiledbsoma/util.format_elapsed"></a>

#### format\_elapsed

```python
def format_elapsed(start_stamp: float, message: str) -> str
```

Returns the message along with an elapsed-time indicator, with end time relative to start
start from `get_start_stamp`. Used for annotating elapsed time of a task.

<a id="tiledbsoma/util.find_csr_chunk_size"></a>

#### find\_csr\_chunk\_size

```python
def find_csr_chunk_size(mat: sp.csr_matrix, start_row_index: int,
                        goal_chunk_nnz: int) -> int
```

Given a CSR matrix and a start row index, returns the number of rows with cumulative NNZ as desired. Context is chunked-COO ingest of larger CSR matrices: if mat is say 8000x9000 but sparse, maybe we'll read rows 0:45 as one chunk and convert that to COO and ingest, then maybe rows 46:78 as a second chunk and convert that to COO and ingest, and so on.

**Arguments**:

- `mat`: The input CSR matrix.
- `start_row_index`: the row index at which to start a chunk.
- `goal_chunk_nnz`: Desired number of non-zero array entries for the chunk.

<a id="tiledbsoma/util.find_csc_chunk_size"></a>

#### find\_csc\_chunk\_size

```python
def find_csc_chunk_size(mat: sp.csc_matrix, start_col_index: int,
                        goal_chunk_nnz: int) -> int
```

Given a CSC matrix and a start column index, returns the number of columns with cumulative nnz as

desired. Context is chunked-COO ingest of larger CSC matrices: if mat is say 8000x9000 but
sparse, maybe we'll read columns 0:45 as one chunk and convert that to COO and ingest, then maybe
columns 46:78 as a second chunk and convert that to COO and ingest, and so on.

**Arguments**:

- `mat`: The input CSC matrix.
- `start_col_index`: the column index at which to start a chunk.
- `goal_chunk_nnz`: Desired number of non-zero array entries for the chunk.

