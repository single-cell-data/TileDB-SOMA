<a id="tiledbsc.raw_group"></a>

# tiledbsc.raw\_group

<a id="tiledbsc.raw_group.RawGroup"></a>

## RawGroup Objects

```python
class RawGroup(TileDBGroup)
```

Nominally for soma raw.

<a id="tiledbsc.raw_group.RawGroup.__init__"></a>

#### \_\_init\_\_

```python
def __init__(uri: str,
             name: str,
             obs: AnnotationDataFrame,
             parent: Optional[TileDBGroup] = None)
```

See the `TileDBObject` constructor.
See `AssayMatrix` for the rationale behind retaining a reference to the `parent_obs` object.

<a id="tiledbsc.raw_group.RawGroup.from_anndata"></a>

#### from\_anndata

```python
def from_anndata(anndata: ad.AnnData) -> None
```

Writes `anndata.raw` to a TileDB group structure.

<a id="tiledbsc.raw_group.RawGroup.to_anndata_raw"></a>

#### to\_anndata\_raw

```python
def to_anndata_raw(obs_labels)
```

Reads TileDB storage and returns the material for an `anndata.Raw` object.
The `obs_labels` must be from the parent object.

