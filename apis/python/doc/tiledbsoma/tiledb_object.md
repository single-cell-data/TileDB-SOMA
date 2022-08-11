<a id="tiledbsc.tiledbsoma/tiledb_object"></a>

# tiledbsc.tiledbsoma/tiledb\_object

<a id="tiledbsc.tiledbsoma/tiledb_object.TileDBObject"></a>

## TileDBObject Objects

```python
class TileDBObject(ABC)
```

Base class for `TileDBArray` and `SOMACollection`.

Manages tiledb_platform_config, context, etc. which are common to both.

<a id="tiledbsc.tiledbsoma/tiledb_object.TileDBObject.__init__"></a>

#### \_\_init\_\_

```python
def __init__(uri: str,
             name: Optional[str] = None,
             *,
             parent: Optional["tiledbsc.tiledbsoma.SOMACollection"] = None,
             tiledb_platform_config: Optional[TileDBPlatformConfig] = None,
             ctx: Optional[tiledb.Ctx] = None)
```

Initialization-handling shared between `TileDBArray` and `SOMACollection`.  Specify
`tiledb_platform_config` and `ctx` for the top-level object; omit them and specify parent for
non-top-level objects. Note that the parent reference is solely for propagating options,
ctx, display depth, etc.

<a id="tiledbsc.tiledbsoma/tiledb_object.TileDBObject.__repr__"></a>

#### \_\_repr\_\_

```python
def __repr__() -> str
```

Fallback string display. Will be overridden by any interesting subclasses.

<a id="tiledbsc.tiledbsoma/tiledb_object.TileDBObject.exists"></a>

#### exists

```python
def exists() -> bool
```

Returns true if the object exists and has the desired class name.

This might be in case an object has not yet been populated, or, if a containing object has
been populated but doesn't have a particular member (e.g. not all `SOMAMeasurement` objects
have a `varp`).

For `tiledb://` URIs this is a REST-server request which we'd like to cache.
However, remove-and-replace use-cases are possible and common in notebooks
and it turns out caching the existence-check isn't a robust approach.

