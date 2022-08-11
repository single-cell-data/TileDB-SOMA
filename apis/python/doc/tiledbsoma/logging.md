<a id="tiledbsc.tiledbsoma/logging"></a>

# tiledbsc.tiledbsoma/logging

<a id="tiledbsc.tiledbsoma/logging.warning"></a>

#### warning

```python
def warning() -> None
```

Sets `tiledbsc.tiledbsoma.logging` to a WARNING level. Use `tiledbsc.tiledbsoma.logging.info()` in notebooks to suppress
progress indicators for data ingestion.

<a id="tiledbsc.tiledbsoma/logging.info"></a>

#### info

```python
def info() -> None
```

Sets `tiledbsc.tiledbsoma.logging` to an INFO level. Use `tiledbsc.tiledbsoma.logging.info()` in notebooks to see
progress indicators for data ingestion.

<a id="tiledbsc.tiledbsoma/logging.debug"></a>

#### debug

```python
def debug() -> None
```

Sets `tiledbsc.tiledbsoma.logging` to an DEBUG level. Use `tiledbsc.tiledbsoma.logging.debug()` in notebooks to see more
detailed progress indicators for data ingestion.

<a id="tiledbsc.tiledbsoma/logging.log_io"></a>

#### log\_io

```python
def log_io(info_message: Optional[str], debug_message: str) -> None
```

Data-ingestion timeframes range widely.  Some folks won't want details for smaller uploads; some
will want details for larger ones.  For I/O and for I/O only, it's helpful to print a short
message at INFO level, or a different, longer message at/beyond DEBUG level.

