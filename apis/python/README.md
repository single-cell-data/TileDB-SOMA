# Overview

This is a Python implementation of the [SOMA API specification](https://github.com/single-cell-data/SOMA/blob/main/abstract_specification.md) for interacting with the [Unified Single-cell Data Model](https://github.com/single-cell-data/SOMA).

# Installation

## Using pip

This code is hosted at [PyPI](https://pypi.org/project/tiledbsoma/), so you can install using Pip:

```shell
$ python -m pip install tiledbsoma
```

To install a specific version:

```shell
$ python -m pip install git+https://github.com/single-cell-data/TileDB-SOMA.git@0.0.6#subdirectory=apis/python
```

To update to the latest version:

```shell
$ python -m pip install --upgrade tiledbsoma
```

## From source

* This requires [`tiledb`](https://github.com/TileDB-Inc/TileDB-Py) (see [./setup.cfg](setup.cfg) for version), in addition to other dependencies in [setup.cfg](./setup.cfg).
* Clone [this repo](https://github.com/single-cell-data/TileDB-SOMA)
* `cd` into your checkout and then `cd apis/python`
* `python -m pip install .`
* Or, if you wish to modify the code and run it, `python -m pip install -v -e .`
* Optionally, if you prefer, you can run that inside `venv`:
  ```shell
  $ python -m venv venv
  $ . ./venv/bin/activate
  $ python -m pip install -v -e .
  ```
* In either case:
  ```shell
  python -m pytest tests
  ```

# Status

Please see [https://github.com/single-cell-data/TileDB-SOMA/issues](https://github.com/single-cell-data/TileDB-SOMA/issues).

# `platform_config` format

When accessing SOMA APIs, TileDB-specific settings can be configured with the [`platform_config`](https://github.com/single-cell-data/SOMA/blob/main/abstract_specification.md#platform-specific-configuration) parameter.
The options accepted by TileDB SOMA are described here, using [TypeScript interface syntax](https://www.typescriptlang.org/docs/handbook/2/objects.html):

```typescript
interface PlatformConfig {
  tiledb?: TDBConfig;
}

interface TDBConfig {
  create?: TDBCreateOptions;
}

interface TDBCreateOptions {
  dims?: { [dim: string]: TDBDimension };
  attrs?: { [attr: string]: TDBAttr };
  allows_duplicates?: bool;

  offsets_filters?: TDBFilter[];
  validity_filters?: TDBFilter[];

  capacity?: number;
  cell_order?: string;
  tile_order?: string;
}

interface TDBDimension {
  filters?: TDBFilter[];
  tile?: number;
}

interface TDBAttr {
  filters?: TDBFilter[];
}

/**
 * Either the name of a filter (in which case it will use
 * the default arguments) or a specification with filter args.
 */
type TDBFilter = string | TDBFilterSpec;

interface TDBFilterSpec {
  /** The name of the filter. */
  _name: string;
  /** kwargs that are passed when constructing the filter. */
  [kwarg: string]: any;
}
```

# Information for developers

Please see the [TileDB-SOMA wiki](https://github.com/single-cell-data/TileDB-SOMA/wiki).
