#!/bin/bash
set -euo pipefail
# Context: .github/workflows/*.yaml
export R_REMOTES_NO_ERRORS_FROM_WARNINGS='true'
Rscript -e 'install.packages(c("remotes", "pkgdown"))'
cd apis/r
Rscript -e 'remotes::install_deps(dependencies = TRUE, upgrade=FALSE)'
# Need dependencies = FALSE to avoid getting latest/sometimes-too-new r-tiledb from CRAN when we
# really want it from our r-universe where we control the version we want:
Rscript -e 'remotes::install_local(dependencies = FALSE)'
cd ../..
