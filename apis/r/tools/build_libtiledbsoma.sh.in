#!/bin/sh

if [ ! -d src/libtiledbsoma ]; then
    echo "No 'src/libtiledbsoma' directory. Exiting."
    exit 1
fi

if [ ! -d src/libtiledbsoma/build-lib ]; then
    mkdir src/libtiledbsoma/build-lib
fi

cwd=`pwd`

cd src/libtiledbsoma/build-lib

## The placeholder is filled in by check_cmake_and_git.R
@cmake@ \
      -DDOWNLOAD_TILEDB_PREBUILT=ON \
      -DTILEDBSOMA_BUILD_CLI=OFF \
      -DTILEDBSOMA_ENABLE_TESTING=OFF \
      -DOVERRIDE_INSTALL_PREFIX=OFF \
      -DCMAKE_INSTALL_PREFIX=${cwd}/inst/tiledbsoma ..

make

make install-libtiledbsoma

rm -rf src/libtiledbsoma/build-lib
