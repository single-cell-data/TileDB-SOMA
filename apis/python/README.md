# Overview

This is a POC Python implementation of the proposed [Unified Single-cell Data Model](https://github.com/single-cell-data/SOMA).

# Installation

## Using pip

This code is hosted at [PyPI](https://pypi.org/project/tiledbsc/), so you can do

```
pip install tiledbsc
```

## From source

* This requires [`tiledb`](https://github.com/TileDB-Inc/TileDB-Py) (see [./setup.cfg](setup.cfg) for version), in addition to other dependencies in [setup.cfg](./setup.cfg).
* Clone [this repo](https://github.com/single-cell-data/TileDB-SingleCell)
* `cd` into your checkout and then `cd apis/python`
* `pip install .` -- or, if you wish to modify the code and run it, `pip install --editable .`
* Optionally, if you prefer, instead:
```
python -m venv venv
. ./venv/bin/activate
pip install .
```
* In either case:

```
python -m pytest tests
```

# Status

Please see [https://github.com/single-cell-data/TileDB-SingleCell/issues](https://github.com/single-cell-data/TileDB-SingleCell/issues).
