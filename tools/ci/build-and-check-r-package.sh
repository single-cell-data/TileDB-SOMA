#!/bin/sh

cd apis/R

echo "::group::Build TileDB R package"
## build package
R CMD build --no-build-vignettes --no-manual .
echo "::endgroup::"

echo "::group::Check TileDB R package"
## set flag to tolerate missing suggested packages (as eg vignette building tools)
export _R_CHECK_FORCE_SUGGESTS_=FALSE
## set an 'under valgrind' variable (see SC-19185)
export _RUNNING_UNDER_VALGRIND_=TRUE
## tell valgrind to use '-s'
export VALGRIND_OPTS="-s"
## check package
R CMD check --use-valgrind --as-cran --no-manual --ignore-vignettes tiledb_*.tar.gz
echo "::endgroup::"

cd ../..
