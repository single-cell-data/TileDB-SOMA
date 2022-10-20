#!/bin/sh

for x in /home/runner/work/TileDB-SOMA/TileDB-SOMA/build/externals/src/ep_tiledb-stamp/ep_tiledb-build-*.log; do
    echo
    echo "================================================================ BEGIN EP BUILD LOG"
    echo $x
    echo "----------------------------------------------------------------"
    cat $x
    echo "================================================================ END   EP BUILD LOG"
done

false # this is invokved as an error-handler so pass error status along
