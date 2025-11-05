# Developer Build Instructions

## System Dependencies

- C++20 compiler
- Python 3.9+
- vcpkg (for dependency management)

Run these commands to setup a fresh Ubuntu 22.04 instance (tested on x86 and Arm):

```
sudo apt update
sudo apt install -y g++ make pkg-config
sudo apt install -y python3 python-is-python3 python3.10-venv python-dev-is-python3
```

______________________________________________________________________

## Clone the source code

```
git clone https://github.com/single-cell-data/TileDB-SOMA.git
cd TileDB-SOMA
```

______________________________________________________________________

## Set up vcpkg

This project uses [vcpkg](https://github.com/microsoft/vcpkg) to manage C++ dependencies (spdlog, fmt, catch2). You have two options:

### Option 1: Use an existing vcpkg installation

Set the `VCPKG_ROOT` environment variable to point to your vcpkg installation:

```bash
export VCPKG_ROOT=/path/to/vcpkg
```

### Option 2: Let CMake fetch vcpkg automatically

If you don't have vcpkg installed, CMake can automatically fetch it for you. Set the `TILEDBSOMA_FETCH_VCPKG` option when configuring:

```bash
cmake -DTILEDBSOMA_FETCH_VCPKG=ON -B build -S libtiledbsoma --preset vcpkg
```

> **Note** - The build script (`scripts/bld`) automatically detects and uses vcpkg if `VCPKG_ROOT` is set in your environment.

______________________________________________________________________

## Set up a Python Virtual Environment

Create a python virtual environment and install the [developer requirements](../apis/python/requirements_dev.txt):

```
python -m venv test/tiledbsoma
source test/tiledbsoma/bin/activate
pip install -r apis/python/requirements_dev.txt
```

______________________________________________________________________

## Python-only Development

Developers who do not need to modify the C++ code must use these build commands:

```
# remove old build artifacts
make clean

# build and install
pip install -v -e apis/python

# test
pytest apis/python
```

This approach leverages the build-system defined in [pyproject.toml](../apis/python/pyproject.toml) to reduce the number of dependencies required to build `tiledbsoma`.

> **Note** - Running `python setup.py develop` is not supported.

______________________________________________________________________

## Python and C++ Development

Developers who plan to modify the C++ code must use these build commands:

```
# clean, build, and install
make install

# test
make test
```

> **Note** - These steps avoid issues when trying to use `cmake` from the [pyproject.toml](../apis/python/pyproject.toml) build-system overlay environment.

> **Note** - All CI and local builds of python, R, and C++ leverage the [../scripts/bld](../scripts/bld) build script, providing a common source for the build flows.

> **Note** - Common developer use cases are captured in the [Makefile](../Makefile) described below.

### Developer Makefile

The [Makefile](../Makefile) automates common developer use cases and promotes sharing of consistent build flows.

```
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

______________________________________________________________________

## Notes

### Dependency Management

This project uses a **hybrid dependency management approach**:

- **vcpkg** provides: `spdlog`, `fmt`, and `catch2` (defined in `vcpkg.json`)
- **Superbuild** provides: `TileDB` (not available in vcpkg, built via ExternalProject)

This hybrid approach ensures:

- Common C++ libraries come from vcpkg for consistency and easier version management
- TileDB is built from source or downloaded as a prebuilt binary via the superbuild system

### System Package Conflicts

If you have system-installed versions of `libfmt`, `libspdlog`, or `tiledb`, they may conflict with the versions required by `tiledbsoma`. The build system prioritizes vcpkg-provided packages, but conflicts can still occur.

**Recommendation**: Uninstall system packages for these dependencies and let the build system manage them:

```bash
# On Ubuntu/Debian
sudo apt remove libfmt-dev libspdlog-dev  # or spdlog, fmt depending on OS
sudo apt remove libtiledb-dev  # if installed
```

As a pro-tip: check the following is gone (and manually remove if necessary):

```bash
/usr/lib/*/cmake/spdlog/spdlogConfig.cmake
```

If you do have reason to have `fmt` and `spdlog` installed system-wide, the following is a known-good configuration on Ubuntu 22.04:

```
$ dpkg -l | egrep "lib(spdlog|fmt)" | cut -c-80
ii  libfmt-dev:amd64                  8.1.1+ds1-2                             am
ii  libfmt8:amd64                     8.1.1+ds1-2                             am
ii  libspdlog-dev:amd64               1:1.9.2+ds-0.2                          am
ii  libspdlog1:amd64                  1:1.9.2+ds-0.2                          am
```

As for the `tiledb` package, if you have reason to have it installed as a separate package, please use the
version matching `libtiledbsoma/cmake/Modules/FindTileDB_EP.cmake`.
