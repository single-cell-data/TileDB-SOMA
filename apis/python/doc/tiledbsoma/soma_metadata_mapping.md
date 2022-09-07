<a id="tiledbsoma/soma_metadata_mapping"></a>

# tiledbsoma/soma\_metadata\_mapping

<a id="tiledbsoma/soma_metadata_mapping.SOMAMetadataMapping"></a>

## SOMAMetadataMapping Objects

```python
class SOMAMetadataMapping()
```

<a id="tiledbsoma/soma_metadata_mapping.SOMAMetadataMapping.get"></a>

#### get

```python
def get(key: str) -> Any
```

Get the value associated with the key.

<a id="tiledbsoma/soma_metadata_mapping.SOMAMetadataMapping.has"></a>

#### has

```python
def has(key: str) -> bool
```

Test for key existence.

<a id="tiledbsoma/soma_metadata_mapping.SOMAMetadataMapping.set"></a>

#### set

```python
def set(key: str, value: Any) -> None
```

Set the value associated with the key.

<a id="tiledbsoma/soma_metadata_mapping.SOMAMetadataMapping.__delete__"></a>

#### \_\_delete\_\_

```python
def __delete__(key: str) -> None
```

Remove the key/value from the collection.

<a id="tiledbsoma/soma_metadata_mapping.SOMAMetadataMapping.__iter__"></a>

#### \_\_iter\_\_

```python
def __iter__() -> Iterator[Any]
```

Iterate over the collection.

<a id="tiledbsoma/soma_metadata_mapping.SOMAMetadataMapping.__len__"></a>

#### \_\_len\_\_

```python
def __len__() -> int
```

Get the length of the map, the number of keys present.

