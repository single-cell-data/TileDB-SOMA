<a id="tiledbsc.soma"></a>

# tiledbsc.soma

<a id="tiledbsc.soma.SOMA"></a>

## SOMA Objects

```python
class SOMA(TileDBGroup)
```

Single-cell group
Class for representing a group of TileDB groups/arrays that constitute an SOMA ('stack of matrices, annotated')
which includes:

* `X` (`AssayMatrixGroup`): a group of one or more labeled 2D sparse arrays that share the same dimensions.
* `obs` (`AnnotationDataframe`): 1D labeled array with column labels for `X`
* `var` (`AnnotationDataframe`): 1D labeled array with row labels for `X`

See also desc-ann.py in this directory for helpful information to
reveal the diversity/variety of HDF5 files we process.

<a id="tiledbsc.soma.SOMA.__init__"></a>

#### \_\_init\_\_

```python
def __init__(uri: str,
             name=None,
             soma_options: Optional[SOMAOptions] = None,
             verbose: Optional[bool] = True,
             config: Optional[tiledb.Config] = None,
             ctx: Optional[tiledb.Ctx] = None,
             parent: Optional[TileDBGroup] = None)
```

Create a new SOMA object. The existing array group is opened at the specified array `uri` if one is present, otherwise a new array group is created.

**Arguments**:

- `uri`: URI of the TileDB group
- `verbose`: Print status messages

<a id="tiledbsc.soma.SOMA.__str__"></a>

#### \_\_str\_\_

```python
def __str__()
```

Implements `print(soma)`.

