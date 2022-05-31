<a id="tiledbsc.soma_collection"></a>

# tiledbsc.soma\_collection

<a id="tiledbsc.soma_collection.SOMACollection"></a>

## SOMACollection Objects

```python
class SOMACollection(TileDBGroup)
```

Implements a collection of `SOMA` objects.

<a id="tiledbsc.soma_collection.SOMACollection.__init__"></a>

#### \_\_init\_\_

```python
def __init__(uri: str,
             name="soco",
             soma_options: Optional[SOMAOptions] = None,
             verbose: Optional[bool] = True,
             config: Optional[tiledb.Config] = None,
             ctx: Optional[tiledb.Ctx] = None,
             parent: Optional[TileDBGroup] = None)
```

Create a new `SOMACollection` object. The existing group is opened at the specified `uri` if one is present, otherwise a new group will be created upon ingest.

**Arguments**:

- `uri`: URI of the TileDB group
- `verbose`: Print status messages

<a id="tiledbsc.soma_collection.SOMACollection.add"></a>

#### add

```python
def add(soma: SOMA) -> None
```

Adds a `SOMA` to the `SOMACollection`.

<a id="tiledbsc.soma_collection.SOMACollection.remove"></a>

#### remove

```python
def remove(soma: SOMA) -> None
```

Removes a `SOMA` from the `SOMACollection`.

<a id="tiledbsc.soma_collection.SOMACollection.__iter__"></a>

#### \_\_iter\_\_

```python
def __iter__() -> List[SOMA]
```

Implements `for soma in soco: ...`

<a id="tiledbsc.soma_collection.SOMACollection.__contains__"></a>

#### \_\_contains\_\_

```python
def __contains__(name: str) -> bool
```

Implements `name in soco`

<a id="tiledbsc.soma_collection.SOMACollection.__getitem__"></a>

#### \_\_getitem\_\_

```python
def __getitem__(name)
```

Returns a `SOMA` element at the given name within the group, or `None` if no such
member exists.  Overloads the `[...]` operator.

