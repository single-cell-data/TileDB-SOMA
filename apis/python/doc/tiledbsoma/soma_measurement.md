<a id="tiledbsc.tiledbsoma/soma_measurement"></a>

# tiledbsc.tiledbsoma/soma\_measurement

<a id="tiledbsc.tiledbsoma/soma_measurement.SOMAMeasurement"></a>

## SOMAMeasurement Objects

```python
class SOMAMeasurement(SOMACollection)
```

A `SOMAMeasurement` is a sub-element of a `SOMAExperiment`, and is otherwise a specialized
`SOMACollection` with pre-defined fields:

`var`: `SOMADataFrame`

Primary annotations on the variable axis, for variables in this measurement (i.e., annotates
columns of `X`). The contents of the `soma_rowid` pseudo-column define the variable index domain,
AKA varid. All variables for this measurement must be defined in this dataframe.

`X`: `SOMACollection` of `SOMASparseNdArray`

A collection of sparse matrices, each containing measured feature values. Each matrix is indexed
by `[obsid, varid]`.

`obsm`: `SOMACollection` of `SOMADenseNdArray`

A collection of dense matrices containing annotations of each `obs` row. Has the same shape as
`obs`, and is indexed with `obsid`.

`obsp`: `SOMACollection` of `SOMASparseNdArray`

A collection of sparse matrices containing pairwise annotations of each `obs` row. Indexed with
`[obsid_1, obsid_2]`.

`varm`: `SOMACollection` of `SOMADenseNdArray`

A collection of dense matrices containing annotations of each `var` row. Has the same shape as
`var`, and is indexed with `varid`.

`varp`: `SOMACollection` of `SOMASparseNdArray`

A collection of sparse matrices containing pairwise annotations of each `var` row. Indexed with
`[varid_1, varid_2]`

<a id="tiledbsc.tiledbsoma/soma_measurement.SOMAMeasurement.__init__"></a>

#### \_\_init\_\_

```python
def __init__(uri: str,
             *,
             name: Optional[str] = None,
             parent: Optional[SOMACollection] = None,
             tiledb_platform_config: Optional[TileDBPlatformConfig] = None,
             ctx: Optional[tiledb.Ctx] = None)
```

Also see the `TileDBObject` constructor.

<a id="tiledbsc.tiledbsoma/soma_measurement.SOMAMeasurement.create"></a>

#### create

```python
def create() -> None
```

Creates the data structure on disk/S3/cloud.

<a id="tiledbsc.tiledbsoma/soma_measurement.SOMAMeasurement.__getattr__"></a>

#### \_\_getattr\_\_

```python
def __getattr__(name: str) -> Any
```

Implements `experiment.var`, `experiment.X`, etc.

<a id="tiledbsc.tiledbsoma/soma_measurement.SOMAMeasurement.constrain"></a>

#### constrain

```python
def constrain() -> None
```

Checks constraints on the collection. Raises an exception if any is violated.

