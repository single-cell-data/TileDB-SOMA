<a id="tiledbsc.v1/util_ann"></a>

# tiledbsc.v1/util\_ann

<a id="tiledbsc.v1/util_ann.describe_ann_file"></a>

#### describe\_ann\_file

```python
def describe_ann_file(input_path: Path,
                      show_summary: bool = True,
                      show_types: bool = False,
                      show_data: bool = False) -> None
```

This is an anndata-describer that goes a bit beyond what `h5ls` does for us.
In particular, it shows us that for one HDF5 file we have `anndata.X` being of type `numpy.ndarray`
while for another HDF5 file we have `anndata.X` being of type `scipy.sparse.csr.csr_matrix`.  This is
crucial information for building I/O logic that accepts a diversity of anndata HDF5 files.

