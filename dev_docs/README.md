# Developer Build Instructions

## System Dependencies

- C++20 compiler
- Python 3.9+ (for Python API)
- R (for R API)
- vcpkg (for dependency management)

## C++ Installation

Both the Python and R API are built on top of an internal C++ `libtiledbsoma` library that is built with CMake. The `libtiledbsoma` library uses [vcpkg](https://github.com/microsoft/vcpkg) to manage C++ dependencies (spdlog, fmt, catch2) except for the core TileDB library which is currently built with a Superbuild. You can also manually install the dependencies using standard CMake search paths, although this method is not actively tested.

_Note on system conflicts: If you have system-installed versions of `libfmt`, `libspdlog`, or `tiledb`, they may conflict with the versions required by `tiledbsoma`. The build system prioritizes vcpkg-provided packages, but conflicts can still occur._

### Dependencies

TileDB core:

You can either manually install TileDB and add it to you CMake module search path (for example, by installing to `/usr/local/`) or you can let `libtiledbsoma` install TileDB-SOMA as part of a Superbuild architecture. If CMake does not find TileDB, it will build the full package itself. Make sure you have an appropriate version of TileDB installed when using a manual install.

vcpkg:

You can either use an existing installation of vcpkg or have CMake download vcpkg as part of the configuration step.

To use an existing vcpkg installation, you can:

1. Set the environment variable `VCPKG_ROOT` to your vcpkg installation. For example, with bash:

   ```
   export VCPKG_ROOT=/path/to/vcpkg
   ```

OR

2. Set the cached CMake variable `CMAKE_TOOLCHAIN_FILE` to the vcpkg toolchain. This can manually be passed to CMake or set with a CMake preset. See the [vcpkg documentation](https://learn.microsoft.com/en-us/vcpkg/users/buildsystems/cmake-integration#cmake_toolchain_file) for more information.

To have CMake install vcpkg as part of the build enable the CMake option with `TILEDBSOMA_FETCH_VCPKG=ON` (this is `OFF` by default). This option will be ignored if either the `CMAKE_TOOLCHAIN_FILE` is set or the `VCPKG_ROOT` environment variable is defined.

### Option 1: Build with CMake (recommended)

The root CMake file for the `libtiledbsoma` library is located in the `libtiledb` subdirectory. This can be called used directly to build and install the C++ library. For example, you can build inside a `build` folder in the root TileDB-SOMA directory using bash by excuting the following commands. _Note: The non-standard re-build and install commands are an artifact of the superbuild architecture which will eventually be replaced._

Configure:

```bash
cmake -S libtiledbsoma -B build [... other options]
```

Build:

```bash
cmake --build build -j $nproc
```

Re-build:

______________________________________________________________________

## Set up a Python Virtual Environment

Create a python virtual environment and install the developer dependencies:

```
python -m venv test/tiledbsoma
source test/tiledbsoma/bin/activate
pip install -e apis/python[dev]
```

Alternatively, you can install the package with all optional dependencies:

```
pip install -e apis/python[all]
```

Test:

```bash
ctest --test-dir build/libtiledbsoma/test
```

Install:

```bash
cmake --build build -t install-libtiledbsoma
```

### Option 2: Build with developer Makefile

The [Makefile](../Makefile) automates common developer use cases. This build process requires an existing installation of vcpkg.

```bash
Usage: make rule [options]

Rules:
  install [options]   Build C++ library and install python module
  r-build [options]   Build C++ static library for R
  update              Incrementally build C++ library and update python module
  test                Run tests
  clean               Remove build artifacts

Options:
  build=BUILD_TYPE    Cmake build type = Release|Debug|RelWithDebInfo|ASAN|TSAN|LSAN|UBSAN|MSAN|Coverage [Release]
  prefix=PREFIX       Install location [dist]
  tiledb=TILEDB_DIST  Absolute path to custom TileDB build

Examples:
  Install Release build

    make install

  Install Debug build of libtiledbsoma and libtiledb

    make install build=Debug

  Install Release build with custom libtiledb

    make install tiledb=$PWD/../TileDB/dist

  Incrementally build C++ changes and update the python module

    make update
```

## Python Installation

The Python developer dependencies are stored in [apis/python/requirements_dev.txt](../apis/python/requirements_dev.txt). These dependencies are necessary for linting and running unit tests. Once the C++ libtiledbsoma library is installed, an editable python build can be installed with pip:

```bash
pip install -e apis/python
```

The test suite can be run with pytest:

```bash
pytest apis/python/tests
```

See the [apis/python/pyproject.toml](../apis/python/pyproject.toml) for available marks for filtering the testing suite.

## R Installation

Once the C++ libtiledbsoma library is installed, an R build can be installed using devtools. From the [apis/r](../apis/r/) directory run the following:

Install dependencies:

```bash
Rscript -e 'remotes::install_deps(".", TRUE)'
```

Initial install:

```bash
Rscript -e 'devtools::install(upgrade=FALSE)'
```

Update:

```bash
Rscript -e 'devtools::load_all(recompile = TRUE)'
```

# Custom vcpkg Triplets for TileDB-SOMA

This directory contains custom vcpkg triplets for building Python wheels with optimized release configurations.

## Purpose

These triplets configure vcpkg to build dependencies (TileDB, spdlog, etc.) in release-only mode with dynamic linking, which is optimal for distributing Python wheels.

## Available Triplets

- `x64-linux-release.cmake` - Linux x86_64 release build
- `arm64-linux-release.cmake` - Linux ARM64 release build
- `x64-osx-release.cmake` - macOS x86_64 (Intel) release build with deployment target 13.0
- `arm64-osx-release.cmake` - macOS ARM64 (Apple Silicon) release build with deployment target 13.0

All triplets are configured with:

- **Build Type**: Release only (no debug symbols, optimized builds)
- **CRT Linkage**: Dynamic (runtime library linked dynamically)
- **Library Linkage**: Dynamic (shared libraries for easier wheel bundling)
- **macOS Deployment Target**: 13.0 (compatibility with macOS Ventura and later)

## Usage

### CI Builds

The triplets are automatically used during CI builds:

1. Included in source distributions (sdist) via `pyproject.toml`
1. Copied to vcpkg during the `before-build` step in `cibuildwheel`
1. Referenced via the `VCPKG_TARGET_TRIPLET` environment variable

### Local Development

For local wheel builds, the triplets are used automatically if building from a source distribution (sdist).

If building from the git repository, you can use them by:

1. Setting `VCPKG_TARGET_TRIPLET` environment variable:

   ```bash
   export VCPKG_TARGET_TRIPLET=x64-linux-release  # or arm64-osx-release, etc.
   ```

1. The build system will automatically copy them to vcpkg during the build

For local development/debugging, you can also use vcpkg's default triplets (e.g., `x64-linux`) which include debug symbols.
