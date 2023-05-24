# Overview

This is the R implementation of the [Unified Single-Cell Data Model](https://github.com/single-cell-data/SOMA).

This `main` branch implements the [updated specfication](https://github.com/single-cell-data/SOMA/blob/main/abstract_specification.md).
Please also see the `main-old` branch for an implementation of the [original specification](https://github.com/single-cell-data/TileDB-SOMA/blob/main-old/spec/specification.md).

# Installation

## Using R-universe

This package is hosted at R-universe so you can do

```shell
> install.packages('tiledbsoma', repos = c('https://tiledb-inc.r-universe.dev',
                                           'https://cloud.r-project.org'))
```

to install it directly (as a pre-made binary on macOS, or from source on Linux) along with all its dependencies.

## From source

### Requirements

* Source installation requires the [`tiledb` R package](https://github.com/TileDB-Inc/TileDB-R) -- which in turn depends on either a local installation of [TileDB Core library](https://github.com/TileDB-Inc/TileDB) or the provided build artifacts.
* In general, source installation of TileDB Core and its packages requires `cmake` and `git` to be installed; these are common tools each operating system provides readily.
* All other R package dependencies are listed in the [DESCRIPTION](https://github.com/single-cell-data/TileDB-SOMA/blob/main/apis/r/DESCRIPTION) file and can be installed via _e.g._
  `remotes::install_deps(".", TRUE)`. In order to build vignettes, `knitr` and `rmarkdown` are required, as is `testthat` for testing. If `testthat` is invoked directly then `pkgbuild` is also needed.
* In addition, this R package also depends on the [`libtiledbsoma` library](https://github.com/single-cell-data/TileDB-SOMA/tree/main/libtiledbsoma) from this repository. It is either installed with the package (as described in the next section), or can be used as a system library (if one is found). In case a system installation is desired, it can provided by following the steps in the [`libtiledbsoma` directory](https://github.com/single-cell-data/TileDB-SOMA/tree/main/libtiledbsoma).


### Step-by-step

* Clone this repository: `git clone https://github.com/single-cell-data/TileDB-SOMA.git`.
* Change into the R API package directory: `cd TileDB-SOMA/apis/r`.
* Optionally, clean the files in the repo: `./cleanup` (this is not needed the first time).
* Optionally, add on optional package with more more test data: `install.packages("pbmc3k.tiledb",  repos="https://ghrr.github.io/drat")`.
* If you have edited any `src/*.cpp` files and changed any function signaturess, then a run `Rscript -e 'Rcpp::compileAttributes()'` will update the `Rcpp`-generated glue code.
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
