# Installation

This requires [`tiledb`](https://github.com/TileDB-Inc/TileDB-Py) (see [./setup.cfg](setup.cfg) for version), in addition to other dependencies in [setup.cfg](./setup.cfg).

After `cd` to `apis/python`:

```
pip install --editable .
```

Optionally, if you prefer, instead:

```
python -m venv venv
. ./venv/bin/activate
pip install .
```

Then:

```
python -m pytest tests
```

# Status

Please see [https://github.com/single-cell-data/TileDB-SingleCell/issues](https://github.com/single-cell-data/TileDB-SingleCell/issues).
