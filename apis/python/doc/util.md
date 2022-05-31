<a id="tiledbsc.util"></a>

# tiledbsc.util

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

