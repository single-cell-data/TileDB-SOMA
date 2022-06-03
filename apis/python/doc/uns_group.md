<a id="tiledbsc.uns_group"></a>

# tiledbsc.uns\_group

<a id="tiledbsc.uns_group.UnsGroup"></a>

## UnsGroup Objects

```python
class UnsGroup(TileDBGroup)
```

Nominally for soma uns.

<a id="tiledbsc.uns_group.UnsGroup.__init__"></a>

#### \_\_init\_\_

```python
def __init__(uri: str, name: str, parent: Optional[TileDBGroup] = None)
```

See the TileDBObject constructor.

<a id="tiledbsc.uns_group.UnsGroup.keys"></a>

#### keys

```python
def keys()
```

For uns, `.keys()` is a keystroke-saver for the more general group-member
accessor `._get_member_names()`.

<a id="tiledbsc.uns_group.UnsGroup.__getitem__"></a>

#### \_\_getitem\_\_

```python
def __getitem__(name)
```

Returns an `UnsArray` or `UnsGroup` element at the given name within the group, or None if
no such member exists.  Overloads the [...] operator.

<a id="tiledbsc.uns_group.UnsGroup.__contains__"></a>

#### \_\_contains\_\_

```python
def __contains__(name)
```

Implements '"namegoeshere" in soma.uns'.

<a id="tiledbsc.uns_group.UnsGroup.__iter__"></a>

#### \_\_iter\_\_

```python
def __iter__() -> List
```

Implements `for element in soma.uns: ...`

<a id="tiledbsc.uns_group.UnsGroup.show"></a>

#### show

```python
def show(display_name="uns")
```

Recursively displays the uns data.

<a id="tiledbsc.uns_group.UnsGroup.from_anndata_uns"></a>

#### from\_anndata\_uns

```python
def from_anndata_uns(uns: ad.compat.OverloadedDict)
```

Populates the uns group for the soma object.

**Arguments**:

- `uns`: anndata.uns.

<a id="tiledbsc.uns_group.UnsGroup.to_dict_of_matrices"></a>

#### to\_dict\_of\_matrices

```python
def to_dict_of_matrices() -> Dict
```

Reads the recursive group/array uns data from TileDB storage and returns them as a recursive dict of matrices.

