# This Makefile captures common developer build and test use cases.

MAKEFLAGS += --no-print-directory

# print help by default
help:

# install
# -------------------------------------------------------------------

# Set default variable values, if non-null
# * build=Debug creates binary artifacts with symbols, e.g. for gdb
# * cmake_verbose=true creates Makefiles that produce full compile lines when executed
build ?= Release
cmake_verbose ?= false

.PHONY: install
install: clean
	@./scripts/bld --prefix=${prefix} --tiledb=${tiledb} --build=${build} --cmake-verbose=${cmake_verbose}
	@TILEDB_PATH=${tiledb} pip install -v -e apis/python

.PHONY: r-build
r-build: clean
	@./scripts/bld --prefix=${prefix} --tiledb=${tiledb} --build=${build} --r-build

# incremental compile and update python install
# -------------------------------------------------------------------
.PHONY: update
update:
	cd build && make -j && make install-libtiledbsoma
	cp dist/lib/lib* apis/python/src/tiledbsoma/

# test
# -------------------------------------------------------------------
.PHONY: test

test: ctest
	pytest apis/python/tests

.PHONY: ctest
ctest: data ctest_update

.PHONY: ctest_update
ctest_update:
	ctest --test-dir build/libtiledbsoma -C Release --verbose --rerun-failed --output-on-failure

.PHONY: data
data:
	@./scripts/prepare-test-data.sh

.PHONY: clean_data
clean_data:
	@./scripts/clean-test-data.sh


# clean
# -------------------------------------------------------------------
.PHONY: clean
clean:
	@rm -rf build dist

.PHONY: cleaner
cleaner:
	@printf "*** dry-run mode: remove -n to actually remove files\n"
	git clean -ffdx -e .vscode -e test/tiledbsoma -n

# help
# -------------------------------------------------------------------
define HELP
Usage: make rule [options]

Rules:
  install [options]   Build C++ library and install python module
  r-build [options]   Build C++ static library with "#define R_BUILD" for R
  update              Incrementally build C++ library and update python module
  test                Run tests
  check-format        Run C++ format check
  format              Run C++ format
  clean               Remove build artifacts

Options:
  build=BUILD_TYPE    Cmake build type = Release|Debug|RelWithDebInfo|ASAN|TSAN|LSAN|UBSAN|MSAN|Coverage [Release]
  prefix=PREFIX       Install location [${PWD}/dist]
  tiledb=TILEDB_DIST  Absolute path to custom TileDB build

Examples:
  Install Release build

    make install

  Install Debug build of libtiledbsoma and libtiledb

    make install build=Debug

  Install Release build with custom libtiledb

    make install tiledb=$$PWD/../TileDB/dist

  Incrementally build C++ changes and update the python module

    make update


endef
export HELP

.PHONY: help
help:
	@printf "$${HELP}"
