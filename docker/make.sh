#!/usr/bin/env bash

set -ex

if ! make install; then
  tdb_log=/TileDB-SOMA/build/externals/src/ep_tiledb-stamp/ep_tiledb-configure-err.log
  if [ -e $tdb_log ]; then
    cat /TileDB-SOMA/build/externals/src/ep_tiledb-stamp/ep_tiledb-configure-err.log
  fi
  vcpkg_log=/TileDB-SOMA/build/externals/src/ep_tiledb-build/vcpkg-bootstrap.log
  if [ -e $vcpkg_log ]; then
    cat $vcpkg_log
  fi
fi
