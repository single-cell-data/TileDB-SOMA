# Overview

This is a R implementation of the [Unified Single-cell Data Model](https://github.com/single-cell-data/SOMA).

This branch, `main`, implements the [updated specfication](https://github.com/single-cell-data/SOMA/blob/main/abstract_specification.md).  Please also see the `main-old` branch which implements the [original specification](https://github.com/single-cell-data/TileDB-SOMA/blob/main-old/spec/specification.md).

# Installation

## Using R-universe

This code is hosted at R-universe, so you can do

```shell
$ install.packages('tiledbsoma', repos = c('https://tiledb-inc.r-universe.dev', 'https://cloud.r-project.org'))
```

## Note on building from source

Note that building from source (using `R CMD INSTALL .`) implies also building
[libtiledbsoma](https://github.com/single-cell-data/TileDB-SOMA/tree/main/libtiledbsoma)
from source so additional requirements listed in its
[README.md](https://github.com/single-cell-data/TileDB-SOMA/blob/main/libtiledbsoma/README.md)
may apply.

# Status

Please see [https://github.com/single-cell-data/TileDB-SOMA/issues](https://github.com/single-cell-data/TileDB-SOMA/issues).

# Information for developers

Please see the [TileDB-SOMA wiki](https://github.com/single-cell-data/TileDB-SOMA/wiki).
