#!/bin/sh

## This script is invoked from the tarball-getting script `get_tarballs.R` which
## itself is called from `configure`.

cwd=`pwd`

## This helps a little
mv inst/tiledbsoma/TileDB-SOMA-* inst/tiledbsoma/TileDB-SOMA

## Now make build/ and call cmake; make; make install
mkdir inst/tiledbsoma/TileDB-SOMA/libtiledbsoma/build && \
    cd inst/tiledbsoma/TileDB-SOMA/libtiledbsoma/build && \
    cmake -DDOWNLOAD_TILEDB_PREBUILT=ON \
          -DTILEDBSOMA_BUILD_CLI=OFF \
          -DTILEDBSOMA_ENABLE_TESTING=OFF \
          -DOVERRIDE_INSTALL_PREFIX=OFF \
          -DCMAKE_INSTALL_PREFIX=${cwd}/inst/tiledbsoma .. && \
    make && \
    make install-libtiledbsoma && \
    cd -

rm -rf inst/tiledbsoma/TileDB-SOMA
