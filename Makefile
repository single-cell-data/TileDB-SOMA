# This Makefile captures common developer build and test use cases.

MAKEFLAGS += --no-print-directory

# print help by default
help:

# install 
# -------------------------------------------------------------------
.PHONY: install
install: clean
	@./scripts/bld --prefix=${prefix} --tiledb=${tiledb} --build=${build}
	@pip install -v -e apis/python

.PHONY: r-build
r-build: clean
	@./scripts/bld --prefix=${prefix} --tiledb=${tiledb} --r-build

# incremental compile and update python install
# -------------------------------------------------------------------
.PHONY: update
update:
	cd build && make -j && make install-libtiledbsoma
	cp dist/lib/* apis/python/src/tiledbsoma/

# test
# -------------------------------------------------------------------
.PHONY: test
test: data
	ctest --test-dir build/libtiledbsoma -C Release --verbose
	pytest

.PHONY: data
data:
	./apis/python/tools/ingestor \
		--ifexists replace \
		--soco \
		-o test/soco \
		-n \
		data/pbmc3k_processed.h5ad \
		data/10x-pbmc-multiome-v1.0/subset_100_100.h5ad

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
  r-build [options]   Build C++ static library for R
  update              Incrementally build C++ library and update python module
  test                Run tests
  clean               Remove build artifacts

Options:
  build=BUILD_TYPE    Cmake build type = Release|Debug|RelWithDebInfo|Coverage [Release]
  prefix=PREFIX       Install location [${PWD}/dist]
  tiledb=TILEDB_DIST  Absolute path to custom TileDB build 

Examples:
  Install Release build

    make install

  Install Debug build of libtiledbsoma and libtiledb

    make install build=Debug

  Install Release build with custom libtiledb

    make install tiledb=$$PWD/../TileDB/dist
    export LD_LIBRARY_PATH=$$PWD/../TileDB/dist/lib:$$LD_LIBRARY_PATH     # linux
    export DYLD_LIBRARY_PATH=$$PWD/../TileDB/dist/lib:$$DYLD_LIBRARY_PATH # macos

  Incrementally build C++ changes and update the python module

    make update


endef 
export HELP

.PHONY: help
help:
	@printf "$${HELP}"

