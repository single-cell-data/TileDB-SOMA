# Overview

This is a R implementation of the [Unified Single-cell Data Model](https://github.com/single-cell-data/SOMA).

This branch, `main`, implements the [updated specfication](https://github.com/single-cell-data/SOMA/blob/main/abstract_specification.md).  Please also see the `main-old` branch which implements the [original specification](https://github.com/single-cell-data/TileDB-SOMA/blob/main-old/spec/specification.md).

# Installation

## Using R-universe

This code is hosted at R-universe, so you can do

```shell
$ install.packages('tiledbsoma', repos = c('https://tiledb-inc.r-universe.dev', 'https://cloud.r-project.org'))
```

## From source

### Requirements

* Source installation requires the [`tiledb` R package](https://github.com/TileDB-Inc/TileDB-R) (which in turn depends on the [`tiledb` Core library](https://github.com/TileDB-Inc/TileDB)).
* This, and the other R package dependencies, are listed in the [DESCRIPTION](https://github.com/single-cell-data/TileDB-SOMA/blob/main/apis/r/DESCRIPTION) file and can be installed via _e.g_ `remotes::install_deps(".")`.
* In addition, the R package also depends on the [`libtiledbsoma` library](https://github.com/single-cell-data/TileDB-SOMA/tree/main/libtiledbsoma) -- which gets installed with the package as
  described in the next section.

### Step-by-step

* Clone this repository: `git clone git@github.com:single-cell-data/TileDB-SOMA.git`
* Change into the R API package directory: `cd TileDB-SOMA/apis/r`
* Optionally, clean the file in to the repo: `./cleanup` (this is not needed the first time)
* Optionally, update the `libtiledbsoma` sources: `./copy_source.sh` (which updates the includes tarball of `libtiledbsoma`).
* Build the R package source tarball from the repository sources: `R CMD build .` (which will also build `libtiledbsoma` from source; other dependencies are required as described in the previous section).
* Install the package from the tarball: `R CMD INSTALL tiledbsoma_*.tar.gz`

Once installed successfully, the package sources can be edited and re-installed.
The optional clean step ensures a full rebuild, and the optional copy of `libtiledbsoma` ensures it is updated too.

# Status

Please see [https://github.com/single-cell-data/TileDB-SOMA/issues](https://github.com/single-cell-data/TileDB-SOMA/issues).

# Information for developers

Please see the [TileDB-SOMA wiki](https://github.com/single-cell-data/TileDB-SOMA/wiki).
