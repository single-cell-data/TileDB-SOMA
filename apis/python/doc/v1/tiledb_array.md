<a id="tiledbsc.v1/tiledb_array"></a>

# tiledbsc.v1/tiledb\_array

<a id="tiledbsc.v1/tiledb_array.TileDBArray"></a>

## TileDBArray Objects

```python
class TileDBArray(TileDBObject)
```

Wraps arrays from TileDB-Py by retaining a URI, options, etc.
Also serves as an abstraction layer to hide TileDB-specific details from the API, unless
requested.

<a id="tiledbsc.v1/tiledb_array.TileDBArray.__init__"></a>

#### \_\_init\_\_

```python
def __init__(uri: str,
             *,
             name: Optional[str] = None,
             parent: Optional["tiledbsc.v1.SOMACollection"] = None)
```

See the TileDBObject constructor.

<a id="tiledbsc.v1/tiledb_array.TileDBArray.get_schema"></a>

#### get\_schema

```python
def get_schema() -> pa.Schema
```

Return data schema, in the form of an Arrow Schema

