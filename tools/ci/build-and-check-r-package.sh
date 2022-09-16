#!/bin/sh

cd apis/R

echo "::group::Build TileDB R package"
## build package
R CMD build --no-build-vignettes --no-manual .
echo "::endgroup::"

echo "::group::Check TileDB R package"
## set flag to tolerate missing suggested packages (as e.g. vignette-building tools)
export _R_CHECK_FORCE_SUGGESTS_=FALSE
## check package
R CMD check --as-cran --no-manual --ignore-vignettes tiledbsoma_*.tar.gz
echo "::endgroup::"

cd ../..
