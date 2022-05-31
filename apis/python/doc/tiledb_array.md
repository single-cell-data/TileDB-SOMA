<a id="tiledbsc.tiledb_array"></a>

# tiledbsc.tiledb\_array

<a id="tiledbsc.tiledb_array.TileDBArray"></a>

## TileDBArray Objects

```python
class TileDBArray(TileDBObject)
```

Wraps arrays from TileDB-Py by retaining a URI, verbose flag, etc.
Also serves as an abstraction layer to hide TileDB-specific details from the API, unless
requested.

<a id="tiledbsc.tiledb_array.TileDBArray.__init__"></a>

#### \_\_init\_\_

```python
def __init__(uri: str, name: str, parent: Optional[TileDBGroup] = None)
```

See the TileDBObject constructor.

<a id="tiledbsc.tiledb_array.TileDBArray.exists"></a>

#### exists

```python
def exists() -> bool
```

Tells whether or not there is storage for the array. This might be in case a SOMA
object has not yet been populated, e.g. before calling `from_anndata` -- or, if the
SOMA has been populated but doesn't have this member (e.g. not all SOMAs have a `varp`).

<a id="tiledbsc.tiledb_array.TileDBArray.tiledb_array_schema"></a>

#### tiledb\_array\_schema

```python
def tiledb_array_schema()
```

Returns the TileDB array schema.

<a id="tiledbsc.tiledb_array.TileDBArray.dim_names"></a>

#### dim\_names

```python
def dim_names() -> List[str]
```

Reads the dimension names from the schema: for example, ['obs_id', 'var_id'].

<a id="tiledbsc.tiledb_array.TileDBArray.dim_names_to_types"></a>

#### dim\_names\_to\_types

```python
def dim_names_to_types() -> Dict[str, str]
```

Returns a dict mapping from dimension name to dimension type.

<a id="tiledbsc.tiledb_array.TileDBArray.attr_names"></a>

#### attr\_names

```python
def attr_names() -> List[str]
```

Reads the attribute names from the schema: for example, the list of column names in a dataframe.

<a id="tiledbsc.tiledb_array.TileDBArray.attr_names_to_types"></a>

#### attr\_names\_to\_types

```python
def attr_names_to_types() -> Dict[str, str]
```

Returns a dict mapping from attribute name to attribute type.

<a id="tiledbsc.tiledb_array.TileDBArray.has_attr_name"></a>

#### has\_attr\_name

```python
def has_attr_name(attr_name: str) -> bool
```

Returns true if the array has the specified attribute name, false otherwise.

<a id="tiledbsc.tiledb_array.TileDBArray.set_soma_object_type_metadata"></a>

#### set\_soma\_object\_type\_metadata

```python
def set_soma_object_type_metadata() -> None
```

This helps nested-structured traversals (especially those that start at the SOMACollection
level) confidently navigate with a minimum of introspection on group contents.

