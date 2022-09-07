<a id="tiledbsoma/io"></a>

# tiledbsoma/io

<a id="tiledbsoma/io.from_h5ad"></a>

#### from\_h5ad

```python
def from_h5ad(experiment: SOMAExperiment, input_path: Path,
              measurement_name: str) -> None
```

Reads an .h5ad file and writes to a TileDB group structure.

<a id="tiledbsoma/io.from_anndata"></a>

#### from\_anndata

```python
def from_anndata(experiment: SOMAExperiment, anndata: ad.AnnData,
                 measurement_name: str) -> None
```

Top-level writer method for creating a TileDB group for a `SOMAExperiment` object.

<a id="tiledbsoma/io.to_h5ad"></a>

#### to\_h5ad

```python
def to_h5ad(experiment: SOMAExperiment, h5ad_path: Path,
            measurement_name: str) -> None
```

Converts the experiment group to anndata format and writes it to the specified .h5ad file.
As of 2022-05-05 this is an incomplete prototype.

<a id="tiledbsoma/io.to_anndata"></a>

#### to\_anndata

```python
def to_anndata(experiment: SOMAExperiment, *,
               measurement_name: str) -> ad.AnnData
```

Converts the experiment group to anndata. Choice of matrix formats is following
what we often see in input .h5ad files:

* X as `scipy.sparse.csr_matrix`
* obs,var as `pandas.dataframe`
* obsm,varm arrays as `numpy.ndarray`
* obsp,varp arrays as `scipy.sparse.csr_matrix`

