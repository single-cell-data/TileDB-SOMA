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
             config: Optional[tiledb.Config] = None,
             ctx: Optional[tiledb.Ctx] = None,
             parent: Optional[TileDBGroup] = None)
```

Create a new `SOMACollection` object. The existing group is opened at the

specified `uri` if one is present, otherwise a new group will be created upon ingest.

**Arguments**:

- `uri`: URI of the TileDB group

<a id="tiledbsc.soma_collection.SOMACollection.__len__"></a>

#### \_\_len\_\_

```python
def __len__() -> int
```

Implements `len(soco)`. Returns the number of elements in the collection.

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

<a id="tiledbsc.soma_collection.SOMACollection.keys"></a>

#### keys

```python
def keys() -> None
```

Returns the names of the SOMAs in the collection.

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
def __getitem__(name) -> SOMA
```

Returns a `SOMA` element at the given name within the group, or `None` if no such
member exists.  Overloads the `[...]` operator.

<a id="tiledbsc.soma_collection.SOMACollection.query"></a>

#### query

```python
def query(obs_attr_names: Optional[List[str]] = None,
          obs_query_string: str = None,
          obs_ids: List[str] = None,
          var_attr_names: Optional[List[str]] = None,
          var_query_string: str = None,
          var_ids: List[str] = None) -> Optional[SOMASlice]
```

Subselects the obs, var, and X/data using the specified queries on obs and var,
concatenating across SOMAs in the collection.  Queries use the TileDB-Py `QueryCondition`
API. If `obs_query_string` is `None`, the `obs` dimension is not filtered and all of `obs`
is used; similiarly for `var`. Return value of `None` indicates an empty slice.
If `obs_ids` or `var_ids` are not `None`, they are effectively ANDed into the query.
For example, you can pass in a known list of `obs_ids`, then use `obs_query_string`
to further restrict the query.

<a id="tiledbsc.soma_collection.SOMACollection.find_unique_obs_values"></a>

#### find\_unique\_obs\_values

```python
def find_unique_obs_values(obs_label: str) -> List
```

Given an `obs` label such as `cell_type` or `tissue`, returns a list of unique values for
that label among all SOMAs in the collection.

<a id="tiledbsc.soma_collection.SOMACollection.find_unique_var_values"></a>

#### find\_unique\_var\_values

```python
def find_unique_var_values(var_label: str) -> List
```

Given an `var` label such as `feature_name`, returns a list of unique values for
that label among all SOMAs in the collection.

