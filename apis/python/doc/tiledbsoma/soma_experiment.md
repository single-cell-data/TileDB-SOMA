<a id="tiledbsoma/soma_experiment"></a>

# tiledbsoma/soma\_experiment

<a id="tiledbsoma/soma_experiment.SOMAExperiment"></a>

## SOMAExperiment Objects

```python
class SOMAExperiment(SOMACollection)
```

`obs`: Primary annotations on the observation axis. The contents of the `soma_rowid` pseudo-column define
the observation index domain, aka `obsid`. All observations for the SOMAExperiment must be
defined in this dataframe.

`ms`: A collection of named measurements.

<a id="tiledbsoma/soma_experiment.SOMAExperiment.__init__"></a>

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

<a id="tiledbsoma/soma_experiment.SOMAExperiment.create"></a>

#### create

```python
def create() -> None
```

Creates the data structure on disk/S3/cloud.

<a id="tiledbsoma/soma_experiment.SOMAExperiment.__getattr__"></a>

#### \_\_getattr\_\_

```python
def __getattr__(name: str) -> Any
```

Implements `experiment.obs` and `experiment.ms`.

<a id="tiledbsoma/soma_experiment.SOMAExperiment.constrain"></a>

#### constrain

```python
def constrain() -> None
```

Checks constraints on the `SOMAExperiment`. Raises an exception if any is violated.

