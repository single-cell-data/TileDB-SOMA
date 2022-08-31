<a id="tiledbsc.soma_slice"></a>

# tiledbsc.soma\_slice

<a id="tiledbsc.soma_slice.SOMASlice"></a>

## SOMASlice Objects

```python
class SOMASlice(TileDBGroup)
```

In-memory-only object for ephemeral extracting out of a SOMA. Can be used to _construct_ a SOMA
but is not a SOMA (which would entail out-of-memory storage).  Nothing more than a collection of
pandas.DataFrame objects. No raw or uns.

<a id="tiledbsc.soma_slice.SOMASlice.__init__"></a>

#### \_\_init\_\_

```python
def __init__(X: Dict[str, Union[pd.DataFrame, pa.Table, Matrix]],
             obs: Union[pd.DataFrame, pa.Table], var: Union[pd.DataFrame,
                                                            pa.Table])
```

Constructs an in-memory `SOMASlice` object. This is a simple collection of obs, var, and X dataframes.

<a id="tiledbsc.soma_slice.SOMASlice.__repr__"></a>

#### \_\_repr\_\_

```python
def __repr__() -> str
```

Default display of SOMASlice.

<a id="tiledbsc.soma_slice.SOMASlice.to_anndata"></a>

#### to\_anndata

```python
def to_anndata() -> ad.AnnData
```

Constructs an `AnnData` object from the current `SOMASlice` object.

<a id="tiledbsc.soma_slice.SOMASlice.concat"></a>

#### concat

```python
@classmethod
def concat(cls, soma_slices: Sequence[SOMASlice]) -> Optional[SOMASlice]
```

Concatenates multiple `SOMASlice` objects into a single one. Implemented using `AnnData`'s
`concat`. Requires that all slices share the same `obs` and `var` keys. Please
see the `SOMA` class method `find_common_obs_and_var_keys`.

