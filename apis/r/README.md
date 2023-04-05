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
* In general, source installation of TileDB and its packages may require `cmake` and `git` to be installed; these are common tools each operating system provides readily.
* Other R package dependencies are listed in the [DESCRIPTION](https://github.com/single-cell-data/TileDB-SOMA/blob/main/apis/r/DESCRIPTION) file and can be installed via _e.g._
  `remotes::install_deps(".", TRUE)`. In order to build vignettes, `knitr` and `rmarkdown` are required as is `testthat` for testing. If `testthat` is invoked directly then `pkgbuild` is also needed.
* In addition, the R package also depends on the [`libtiledbsoma` library](https://github.com/single-cell-data/TileDB-SOMA/tree/main/libtiledbsoma) -- which gets installed with the package as
  described in the next section.

### Step-by-step

* Clone this repository: `git clone https://github.com/single-cell-data/TileDB-SOMA.git`
* Change into the R API package directory: `cd TileDB-SOMA/apis/r`
* Optionally, clean the files in the repo: `./cleanup` (this is not needed the first time)
* Access to more test data: `install.packages("pbmc3k.tiledb",  repos="https://ghrr.github.io/drat")`
`
* Optionally, update the `libtiledbsoma` sources: `./copy_source.sh` (which updates the includes tarball of `libtiledbsoma`).
* If you have edited any `src/*.cpp` with `RcppExport` then run `Rscript -e 'Rcpp::compileAttributes()'`
* Build the R package source tarball from the repository sources: `R CMD build .` (which will also build `libtiledbsoma` from source; other dependencies are required as described in the previous section)
* If you have changed any signatures, run `Rscript -e 'roxygen2::roxygenise()'`
* Check and test the package from the tarball: `R CMD check --no-vignettes --no-manual tiledbsoma_*.tar.gz`
  * For quicker iteration, run `Rscript -e 'testthat::test_local("tests/testthat")'`
* Install the package from the tarball: `R CMD INSTALL tiledbsoma_*.tar.gz`

Once installed successfully, the package sources can be edited and re-installed.
The optional clean step ensures a full rebuild, and the optional copy of `libtiledbsoma` ensures it is updated too.

# Status

Please see [https://github.com/single-cell-data/TileDB-SOMA/issues](https://github.com/single-cell-data/TileDB-SOMA/issues).

# Information for developers

Please see the [TileDB-SOMA wiki](https://github.com/single-cell-data/TileDB-SOMA/wiki).
