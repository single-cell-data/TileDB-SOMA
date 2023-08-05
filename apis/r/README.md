# Overview

This is the R implementation of the [SOMA API specification](https://github.com/single-cell-data/SOMA/blob/main/abstract_specification.md).

# Installation

## Release packages

TileDB-SOMA releases are available on R-universe and [Conda](https://anaconda.org/tiledb/r-tiledbsoma), and can be installed directly from R or `mamba` as indicated below.

```r
install.packages('tiledbsoma', repos = c('https://tiledb-inc.r-universe.dev',
                                         'https://cloud.r-project.org'))
```

System prerequisites include `cmake` and `git`.

```bash
mamba install -c conda-forge -c tiledb r-tiledbsoma
```

The r-universe repo serves macOS binaries and the source package for other Unix-like platforms. The conda channel serves binaries for multiple architectures.

## From source

To install the very latest tiledbsoma development version (our `main` branch), use [`remotes::install_github()`](https://cran.r-project.org/web/packages/remotes/readme/README.html):

```r
remotes::install_github("https://github.com/single-cell-data/TileDB-SOMA", subdir = "apis/r")
```

### Requirements

* Source installation requires the [`tiledb` R package](https://github.com/TileDB-Inc/TileDB-R) -- which in turn depends on either a local installation of [TileDB Core library](https://github.com/TileDB-Inc/TileDB) or the provided build artifacts.
* In general, source installation of TileDB Core and its packages requires `cmake` and `git` to be installed; these are common tools each operating system provides readily.
* All other R package dependencies are listed in the [DESCRIPTION](https://github.com/single-cell-data/TileDB-SOMA/blob/main/apis/r/DESCRIPTION) file and can be installed via _e.g._
  `remotes::install_deps(".", dependencies=TRUE)` which will also install suggested packages. In order to build vignettes, `knitr` and `rmarkdown` are required (and will be installed), as is `testthat` for testing. If `testthat` is invoked directly then `pkgbuild` is also needed (but is not installed as not listed in `DESCRIPTION`).
* In addition, this R package also depends on the [`libtiledbsoma` library](https://github.com/single-cell-data/TileDB-SOMA/tree/main/libtiledbsoma) from this repository. It is either installed with the package (as described in the next section), or can be used as a system library (if one is found). A system installation can be provided by following the steps in the [`libtiledbsoma` directory](https://github.com/single-cell-data/TileDB-SOMA/tree/main/libtiledbsoma).


### Development setup

To set up and install from a local clone of this git repository:

* Clone this repository: `git clone https://github.com/single-cell-data/TileDB-SOMA.git`.
* Change into the R API package directory: `cd TileDB-SOMA/apis/r`.
* Optionally, clean the files in the repo: `./cleanup` (this is not needed the first time).
* Optionally, add on optional package with more more test data: `install.packages("pbmc3k.tiledb",  repos="https://ghrr.github.io/drat")`.
* If you have edited any `src/*.cpp` files and changed any function signatures, then running `Rscript -e 'Rcpp::compileAttributes()'` will update the `Rcpp`-generated glue code.
* If you have changed any C++ function header documentation, run `Rscript -e 'roxygen2::roxygenise()'` to update the corresponding R files.
* Build the R package source tarball from the repository sources: `R CMD build .` (which will also include the `libtiledbsoma` source via a repository soft-link); other dependencies are required as described in the previous section).
* Optionally, check and test the package from the tarball skipping vignettes and manuals (which need `texlive` or equivalent): `R CMD check --no-vignettes --no-manual tiledbsoma_*.tar.gz`.
* Finally, install the package from the tarball: `R CMD INSTALL tiledbsoma_*.tar.gz`.  During this installation presence of the two C++ libraries (TileDB Core, TileDB-SOMA) is tested for. TileDB Core builds can be downloaded as needed, and the TileDB-SOMA library is built if needed. (We plan to provide a downloadable artifact for it too.)

Once installed successfully, the package sources can be edited and re-installed iteratively.
The optional clean step ensures a full rebuild, and the optional copy of `libtiledbsoma` ensures it is updated too.

# Status

Please see [https://github.com/single-cell-data/TileDB-SOMA/issues](https://github.com/single-cell-data/TileDB-SOMA/issues).

# Information for developers

Please see the [TileDB-SOMA wiki](https://github.com/single-cell-data/TileDB-SOMA/wiki).

This `main` branch implements the [updated specification](https://github.com/single-cell-data/SOMA/blob/main/abstract_specification.md). Please also see the `main-old` branch for an implementation of the [original specification](https://github.com/single-cell-data/TileDB-SOMA/blob/main-old/spec/specification.md).
