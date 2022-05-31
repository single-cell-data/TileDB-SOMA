<a id="tiledbsc.tiledb_object"></a>

# tiledbsc.tiledb\_object

<a id="tiledbsc.tiledb_object.TileDBObject"></a>

## TileDBObject Objects

```python
class TileDBObject()
```

Base class for `TileDBArray` and `TileDBGroup`. Manages soma_options, context, etc. which are common
to both.

<a id="tiledbsc.tiledb_object.TileDBObject.__init__"></a>

#### \_\_init\_\_

```python
def __init__(uri: str,
             name: str,
             parent=None,
             soma_options: Optional[SOMAOptions] = None,
             verbose: Optional[bool] = True,
             ctx: Optional[tiledb.Ctx] = None)
```

Initialization-handling shared between `TileDBArray` and `TileDBGroup`.  Specify soma_options,
verbose, and ctx for the top-level object; omit them and specify parent for non-top-level
objects. Note that the parent reference is solely for propagating options, ctx, display
depth, etc.

