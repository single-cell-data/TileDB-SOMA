# This Makefile captures common developer build and test use cases.

MAKEFLAGS += --no-print-directory

# - incrementally compile and install the c++ code
# - copy the shared objects on top of the python dev build
.PHONY: update
update:
	@cd build && make -j && make install-libtiledbsoma
	cp dist/lib/* apis/python/src/tiledbsoma/

# - run the clean rule
# - build the c++ code using the environment's cmake (not the python build-system cmake)
#   - builds a Release with Debug Info build
# - install the python in editable mode (dev build)
.PHONY: install
install: clean
	./scripts/bld
	pip install -v -e apis/python

# - run the clean rule
# - build the c++ code using the environment's cmake (not the python build-system cmake)
#   - builds a Debug build
# - install the python in editable mode (dev build)
.PHONY: install-debug
install-debug: clean
	BUILD_TYPE=Debug ./scripts/bld
	pip install -v -e apis/python

# - run the clean rule
# - build the c++ static library for r
.PHONY: r-build
r-build: clean
	R_BUILD=1 ./scripts/bld

# - run c++ unit tests
# - run pytests
.PHONY: test
test: test/soco
	ctest --test-dir build/libtiledbsoma -C Release --verbose
	pytest

# install the test data
test/soco:
	./apis/python/tools/ingestor \
		--ifexists replace \
		--soco \
		-o test/soco \
		-n \
		data/pbmc3k_processed.h5ad \
		data/10x-pbmc-multiome-v1.0/subset_100_100.h5ad

# - remove the c++ build and dist, and the test data
.PHONY: clean
clean:
	rm -rf build dist test/soco

# - run git clean (in dry-run mode)
.PHONY: cleaner
cleaner:
	@echo "*** running in dry-run mode ***"
	git clean -nffdx -e .vscode -e test/tiledbsoma
