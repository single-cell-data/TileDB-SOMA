<a id="tiledbsc.tiledb_group"></a>

# tiledbsc.tiledb\_group

<a id="tiledbsc.tiledb_group.TileDBGroup"></a>

## TileDBGroup Objects

```python
class TileDBGroup(TileDBObject)
```

Wraps groups from TileDB-Py by retaining a URI, options, etc.

<a id="tiledbsc.tiledb_group.TileDBGroup.__init__"></a>

#### \_\_init\_\_

```python
def __init__(uri: str,
             name: str,
             *,
             parent=None,
             soma_options: Optional[SOMAOptions] = None,
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

<a id="tiledbsc.tiledb_group.TileDBGroup.create_unless_exists"></a>

#### create\_unless\_exists

```python
def create_unless_exists() -> None
```

Creates the TileDB group data structure on disk/S3/cloud, unless it already exists.

<a id="tiledbsc.tiledb_group.TileDBGroup.get_object_type"></a>

#### get\_object\_type

```python
def get_object_type() -> str
```

Returns the class name associated with the group.

<a id="tiledbsc.tiledb_group.TileDBGroup.show_metadata"></a>

#### show\_metadata

```python
def show_metadata(recursively=True, indent="") -> None
```

Shows metadata for the group, recursively by default.

