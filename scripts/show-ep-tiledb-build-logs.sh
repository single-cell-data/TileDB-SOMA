#!/bin/sh

for x in /home/runner/work/TileDB-SOMA/TileDB-SOMA/build/externals/src/ep_tiledb-stamp/ep_tiledb-build-*.log; do
    echo
    echo ------------------ $x
    cat $x
done
