<a id="tiledbsc.v1/soma_collection"></a>

# tiledbsc.v1/soma\_collection

<a id="tiledbsc.v1/soma_collection.SOMACollection"></a>

## SOMACollection Objects

```python
class SOMACollection(TileDBObject)
```

Contains a key-value mapping where the keys are string names and the values are any SOMA-defined
foundational or composed type, including SOMACollection, SOMADataFrame, SOMADenseNdArray,
SOMASparseNdArray or SOMAExperiment.

<a id="tiledbsc.v1/soma_collection.SOMACollection.__init__"></a>

#### \_\_init\_\_

```python
def __init__(uri: str,
             *,
             name: Optional[str] = None,
             parent: Optional[SOMACollection] = None,
             tiledb_platform_config: Optional[TileDBPlatformConfig] = None,
             ctx: Optional[tiledb.Ctx] = None)
```

Also see the :class:`TileDBObject` constructor.

<a id="tiledbsc.v1/soma_collection.SOMACollection.create"></a>

#### create

```python
def create() -> None
```

Creates the data structure on disk/S3/cloud.

<a id="tiledbsc.v1/soma_collection.SOMACollection.__len__"></a>

#### \_\_len\_\_

```python
def __len__() -> int
```

Returns the number of members in the collection.  Implements Python's `len(collection)`.

<a id="tiledbsc.v1/soma_collection.SOMACollection.__contains__"></a>

#### \_\_contains\_\_

```python
def __contains__(member_name: str) -> bool
```

Tests for the existence of key in collection.
Implements the `in` operator.

<a id="tiledbsc.v1/soma_collection.SOMACollection.get"></a>

#### get

```python
def get(member_name: str) -> TileDBObject
```

Get the member object associated with the key

<a id="tiledbsc.v1/soma_collection.SOMACollection.keys"></a>

#### keys

```python
def keys() -> Sequence[str]
```

Gets the names of the members of the collection.

<a id="tiledbsc.v1/soma_collection.SOMACollection.__getitem__"></a>

#### \_\_getitem\_\_

```python
def __getitem__(member_name: str) -> TileDBObject
```

Get the member object associated with the key, when invoked as `collection["namegoeshere"]`.

<a id="tiledbsc.v1/soma_collection.SOMACollection.__getattr__"></a>

#### \_\_getattr\_\_

```python
def __getattr__(member_name: str) -> TileDBObject
```

Get the member object associated with the key, when invoked as `collection.namegoeshere`.

<a id="tiledbsc.v1/soma_collection.SOMACollection.set"></a>

#### set

```python
def set(member: TileDBObject, *, relative: Optional[bool] = None) -> None
```

Adds a member to the collection.

<a id="tiledbsc.v1/soma_collection.SOMACollection.delete"></a>

#### delete

```python
def delete(member_name: str) -> None
```

Removes a member from the collection, when invoked as `collection.delete("namegoeshere")`.

<a id="tiledbsc.v1/soma_collection.SOMACollection.__delattr__"></a>

#### \_\_delattr\_\_

```python
def __delattr__(member_name: str) -> None
```

Removes a member from the collection, when invoked as `del collection.namegoeshere`.

<a id="tiledbsc.v1/soma_collection.SOMACollection.__delitem__"></a>

#### \_\_delitem\_\_

```python
def __delitem__(member_name: str) -> None
```

Removes a member from the collection, when invoked as `del collection["namegoeshere"]`.

<a id="tiledbsc.v1/soma_collection.SOMACollection.__iter__"></a>

#### \_\_iter\_\_

```python
def __iter__() -> Iterator[Any]
```

Iterates over the collection.  Implements Python `for member in collection: ...` syntax.

<a id="tiledbsc.v1/soma_collection.SOMACollection.__repr__"></a>

#### \_\_repr\_\_

```python
def __repr__() -> str
```

Default display for `SOMACollection`.

