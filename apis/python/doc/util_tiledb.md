<a id="tiledbsc.util_tiledb"></a>

# tiledbsc.util\_tiledb

<a id="tiledbsc.util_tiledb.show_soma_schemas"></a>

#### show\_soma\_schemas

```python
def show_soma_schemas(soma_uri: str, ctx: Optional[tiledb.Ctx] = None) -> None
```

Show some summary information about an ingested TileDB Single-Cell Group.  This tool goes a bit beyond `print(tiledb.group.Group(soma_uri))` by also revealing array schema. Additionally, by employing encoded domain-specific knowleldge, it traverses items in the familiar order `X`, `obs`, `var`, etc. rather than using the general-purpose tiledb-group-display function.

<a id="tiledbsc.util_tiledb.show_tiledb_group_array_schemas"></a>

#### show\_tiledb\_group\_array\_schemas

```python
def show_tiledb_group_array_schemas(uri: str,
                                    ctx: Optional[tiledb.Ctx] = None) -> None
```

Recursively show array schemas within a TileDB Group. This function is not specific to
single-cell matrix-API data, and won't necessarily traverse items in a familiar
application-specific order.

