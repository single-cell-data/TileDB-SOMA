# tiledbsoma (R)

R implementation of the [SOMA API specification](https://github.com/single-cell-data/SOMA/blob/main/abstract_specification.md).

## Installation

### From R-universe (recommended)

R-universe serves macOS binaries and source packages for other Unix-like platforms.

```r
install.packages("tiledbsoma", repos = c("https://tiledb-inc.r-universe.dev",
                                         "https://cloud.r-project.org"))
```

Installing from source on Unix-like platforms requires `cmake`, `git`, and a C++ compiler with C++20 support (e.g., `g++` 13+).

### From Conda

Binary packages for multiple architectures are available from [Conda](https://anaconda.org/tiledb/r-tiledbsoma):

```bash
conda install -c conda-forge -c tiledb r-tiledbsoma
```

### From GitHub (development)

To install the latest development version from the `main` branch:

```r
remotes::install_github("https://github.com/single-cell-data/TileDB-SOMA", subdir = "apis/r")
```

### Requirements

Source installation requires:

- [`tiledb` R package](https://github.com/TileDB-Inc/TileDB-R), which depends on a local installation of [TileDB Core](https://github.com/TileDB-Inc/TileDB) or its provided build artifacts
- `cmake` and `git`
- [`libtiledbsoma`](https://github.com/single-cell-data/TileDB-SOMA/tree/main/libtiledbsoma), which is built automatically during package installation or can be provided as a system library (see the [`libtiledbsoma` directory](https://github.com/single-cell-data/TileDB-SOMA/tree/main/libtiledbsoma) for manual installation)

All other R package dependencies are listed in [DESCRIPTION](https://github.com/single-cell-data/TileDB-SOMA/blob/main/apis/r/DESCRIPTION) and can be installed with:

```r
remotes::install_deps(".", dependencies = TRUE)
```

### Troubleshooting: `fmt` and `spdlog` conflicts

If you have system versions of `fmt` or `spdlog` installed, they can conflict with the versions vendored by TileDB.

- **`undefined symbol: _ZTIN3fmt2v912format_errorE`**: Uninstall your system versions of `fmt` and `spdlog`. On Linux, verify removal with `dpkg -l | egrep "lib(spdlog|fmt)"`.
- **`fatal error: spdlog/spdlog.h: No such file or directory`**: Remove the leftover cmake config directory at `/usr/local/lib/cmake/spdlog`, which the system uninstall of `spdlog` fails to clean up.

## Development

### Setup

1. Clone the repository:

   ```bash
   git clone https://github.com/single-cell-data/TileDB-SOMA.git
   cd TileDB-SOMA/apis/r
   ```

1. Install R package dependencies:

   ```r
   remotes::install_deps(".", dependencies = TRUE)
   ```

1. Optionally, install additional data packages used by tests:

   - **`pbmc3k`** provides single-cell data for SingleCellExperiment and SummarizedExperiment ingestion tests:

     ```r
     install.packages("pbmc3k", repos = "https://tiledb-inc.r-universe.dev")
     ```

   - **`pbmc3k.tiledb`** provides pre-built SOMA test data:

     ```r
     install.packages("pbmc3k.tiledb", repos = "https://ghrr.github.io/drat")
     ```

### Build and install

```bash
# Clean previous build artifacts (not needed on first build)
./cleanup

# Build the source tarball (includes libtiledbsoma via symlink)
R CMD build .

# Check the package (skipping vignettes and manual, which require texlive)
R CMD check --no-vignettes --no-manual tiledbsoma_*.tar.gz

# Install
R CMD INSTALL tiledbsoma_*.tar.gz
```

During installation, TileDB Core is downloaded if needed and `libtiledbsoma` is built if a system installation is not found.

### After editing source files

- **Changed `src/*.cpp` function signatures**: Run `Rscript -e 'Rcpp::compileAttributes()'` to update Rcpp-generated glue code.
- **Changed C++ function header documentation**: Run `Rscript -e 'roxygen2::roxygenise()'` to update the corresponding R files.

Once installed, the package sources can be edited and re-installed iteratively. Running `./cleanup` before rebuilding ensures a full rebuild including `libtiledbsoma`.
