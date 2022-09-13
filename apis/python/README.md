# Overview

This is a POC Python implementation of the proposed [Unified Single-cell Data Model](https://github.com/single-cell-data/SOMA).

This branch conforms to [this version](https://github.com/single-cell-data/TileDB-SOMA/blob/main/spec/specification.md) of the specification -- for the [latest version](https://github.com/single-cell-data/SOMA/blob/spec-revision/brainstorming.md) please see the `main` branch of this repository.

# Installation

## Using pip

This code is hosted at [PyPI](https://pypi.org/project/tiledbsc/), so you can do

```
python -m pip install tiledbsc
```

To install a specific version:

```
python -m pip install git+https://github.com/single-cell-data/TileDB-SOMA.git@0.0.6#subdirectory=apis/python
```

To update to the latest version:

```
python -m pip install --upgrade tiledbsc
```

## From source

* This requires [`tiledb`](https://github.com/TileDB-Inc/TileDB-Py) (see [./setup.cfg](setup.cfg) for version), in addition to other dependencies in [setup.cfg](./setup.cfg).
* Clone [this repo](https://github.com/single-cell-data/TileDB-SOMA)
* `cd` into your checkout and then `cd apis/python`
* `python -m pip install .`
* Or, if you wish to modify the code and run it, `python setup.py develop`
* Optionally, if you prefer, you can run that inside `venv`:
```
python -m venv venv
. ./venv/bin/activate
python -m pip install .
```
* In either case:

```
python -m pytest tests
```

# Status

Please see [https://github.com/single-cell-data/TileDB-SOMA/issues](https://github.com/single-cell-data/TileDB-SOMA/issues).
