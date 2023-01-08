#!/bin/bash
set -euo pipefail
# Context: .github/workflows/*.yaml
export R_REMOTES_NO_ERRORS_FROM_WARNINGS='true'
Rscript -e 'install.packages("remotes")'
cd apis/r
Rscript -e 'remotes::install_deps(dependencies = TRUE)'
cd ../..
