# Overview

This is a POC Python implementation of the proposed [Unified Single-cell Data Model](https://github.com/single-cell-data/SOMA).

# Installation

## Using pip

This code is hosted at [PyPI](https://pypi.org/project/tiledbsc/), so you can do

```
python -m pip install tiledbsc
```

To install a specific version:

```
python -m pip install git+https://github.com/single-cell-data/TileDB-SingleCell.git@0.0.6#subdirectory=apis/python
```

To update to the latest version:

```
python -m pip install --upgrade tiledbsc
```

## From source

* This requires [`tiledb`](https://github.com/TileDB-Inc/TileDB-Py) (see [./setup.cfg](setup.cfg) for version), in addition to other dependencies in [setup.cfg](./setup.cfg).
* Clone [this repo](https://github.com/single-cell-data/TileDB-SingleCell)
* `cd` into your checkout and then `cd apis/python`
* `python -m pip install .`
* Optionally, if you prefer `venv`:
```
python -m venv venv
. ./venv/bin/activate
python -m pip install .
```
* If you wish to modify the code and run it in place without reinstalls (this has been done successfully using Python 3.9 and Pip 22.1):
```
python -m venv venv
. ./venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e .
```

* In either case:

```
python -m pytest tests
```

# Status

Please see [https://github.com/single-cell-data/TileDB-SingleCell/issues](https://github.com/single-cell-data/TileDB-SingleCell/issues).
