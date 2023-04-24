#!/bin/sh

## part one: if we are inside the git repo, simply copy over
test -d ../../libtiledbsoma && cp -a ../../libtiledbsoma src/ && cp -a ../../scripts/bld src/libtiledbsoma && rm -rf src/libtiledbsoma/build src/libtiledbsoma/test/__pycache__

## part two: also create a tar.gz 'for keeps' as we may need to install from partial source tree
#test -d ../../libtiledbsoma && cd ../.. && tar -cz -f apis/r/libtiledbsoma.tar.gz libtiledbsoma
