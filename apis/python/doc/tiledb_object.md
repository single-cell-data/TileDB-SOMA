<a id="tiledbsc.tiledb_object"></a>

# tiledbsc.tiledb\_object

<a id="tiledbsc.tiledb_object.TileDBObject"></a>

## TileDBObject Objects

```python
class TileDBObject(ABC)
```

Base class for `TileDBArray` and `TileDBGroup`.

Manages soma_options, context, etc. which are common to both.

<a id="tiledbsc.tiledb_object.TileDBObject.__init__"></a>

#### \_\_init\_\_

```python
def __init__(uri: str,
             name: str,
             *,
             parent: Optional["tiledbsc.TileDBGroup"] = None,
             soma_options: Optional[SOMAOptions] = None,
             ctx: Optional[tiledb.Ctx] = None)
```

Initialization-handling shared between `TileDBArray` and `TileDBGroup`.  Specify soma_options
and ctx for the top-level object; omit them and specify parent for non-top-level
objects. Note that the parent reference is solely for propagating options, ctx, display
depth, etc.

<a id="tiledbsc.tiledb_object.TileDBObject.metadata"></a>

#### metadata

```python
def metadata() -> Mapping[str, Any]
```

Returns metadata from the group/array as a dict.

<a id="tiledbsc.tiledb_object.TileDBObject.has_metadata"></a>

#### has\_metadata

```python
def has_metadata(key: str) -> bool
```

Returns whether metadata is associated with the group/array.

<a id="tiledbsc.tiledb_object.TileDBObject.metadata_keys"></a>

#### metadata\_keys

```python
def metadata_keys() -> Sequence[str]
```

Returns metadata keys associated with the group/array.

<a id="tiledbsc.tiledb_object.TileDBObject.get_metadata"></a>

#### get\_metadata

```python
def get_metadata(key: str) -> Any
```

Returns metadata associated with the group/array.
Raises `KeyError` if there is no such key in the metadata.

<a id="tiledbsc.tiledb_object.TileDBObject.set_metadata"></a>

#### set\_metadata

```python
def set_metadata(key: str, value: Any) -> None
```

Returns metadata associated with the group/array.

<a id="tiledbsc.tiledb_object.TileDBObject.get_object_type"></a>

#### get\_object\_type

```python
def get_object_type() -> str
```

Returns the class name associated with the group/array.

