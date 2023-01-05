#!/bin/bash
set -euo pipefail
# Context: .github/workflows/*.yaml
export R_REMOTES_NO_ERRORS_FROM_WARNINGS='true'

# TODO: Remove this once we have either a TileDB-R 0.17.1 depending on core 2.13.1, or, a TileDB-R
# 0.18 depending on core 2.14.
#
# See also https://github.com/single-cell-data/TileDB-SOMA/issues/654
#
# * Enable repository from tiledb-inc
# * Download and install tiledb in R
Rscript -e '
    options(repos = c(tiledbinc = "https://tiledb-inc.r-universe.dev", CRAN = "https://cloud.r-project.org"));
    install.packages("tiledb")
'

Rscript -e 'install.packages("remotes")'

cd apis/r
Rscript -e 'remotes::install_deps(dependencies = TRUE)'
cd ../..
