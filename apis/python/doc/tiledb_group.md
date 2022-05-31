<a id="tiledbsc.tiledb_group"></a>

# tiledbsc.tiledb\_group

<a id="tiledbsc.tiledb_group.TileDBGroup"></a>

## TileDBGroup Objects

```python
class TileDBGroup(TileDBObject)
```

Wraps groups from TileDB-Py by retaining a URI, verbose flag, etc.

<a id="tiledbsc.tiledb_group.TileDBGroup.__init__"></a>

#### \_\_init\_\_

```python
def __init__(uri: str,
             name: str,
             parent=None,
             soma_options: Optional[SOMAOptions] = None,
             verbose: Optional[bool] = True,
             ctx: Optional[tiledb.Ctx] = None)
```

See the TileDBObject constructor.

<a id="tiledbsc.tiledb_group.TileDBGroup.exists"></a>

#### exists

```python
def exists() -> bool
```

Tells whether or not there is storage for the group. This might be in case a SOMA
object has not yet been populated, e.g. before calling `from_anndata` -- or, if the
SOMA has been populated but doesn't have this member (e.g. not all SOMAs have a `varp`).

