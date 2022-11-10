# Developer Build Instructions

## System Dependencies

* C++17 compiler
* Python 3.8+

Run these commands to setup a fresh Ubuntu 22.04 instance (tested on x86 and Arm):
```
sudo apt update
sudo apt install -y g++ make pkg-config
sudo apt install -y python3 python-is-python3 python3.10-venv python-dev-is-python3
```
---
## Clone the source code
```
git clone https://github.com/single-cell-data/TileDB-SOMA.git
cd TileDB-SOMA
```
---
## Set up a Python Virtual Environment
Create a python virtual environment and install the [developer requirements](../apis/python/requirements_dev.txt):
```
python -m venv test/tiledbsoma
source test/tiledbsoma/bin/activate
pip install -r apis/python/requirements_dev.txt
```
---
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

---

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
  r-build [options]   Build C++ static library with "#define R_BUILD" for R
  update              Incrementally build C++ library and update python module
  test                Run tests
  clean               Remove build artifacts

Options:
  build=BUILD_TYPE    Cmake build type = Release|Debug|RelWithDebInfo|Coverage [Release]
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

---

## Notes

A build system with `libfmt` and `libspdlog` installed may conflict with the required versions for `tiledbsoma`. If that is the case, this configuration is known-good on Ubuntu 22.04:

```
sudo apt install libfmt-dev libspdlog-dev
```

Check:

```
$ dpkg -l | egrep "lib(spdlog|fmt)" | cut -c-80
ii  libfmt-dev:amd64                  8.1.1+ds1-2                             am
ii  libfmt8:amd64                     8.1.1+ds1-2                             am
ii  libspdlog-dev:amd64               1:1.9.2+ds-0.2                          am
ii  libspdlog1:amd64                  1:1.9.2+ds-0.2                          am
```

Similarly, for MacOS:

* `brew uninstall spdlog` (if you have it)
* `brew uninstall fmt` (if you have it)
* `brew uninstall tiledb` (if you have it from `brew`)

All these will be installed for you using compatible versions.
