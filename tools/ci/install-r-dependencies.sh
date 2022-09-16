#!/bin/sh
# Context: .github/workflows/*.yaml
export R_REMOTES_NO_ERRORS_FROM_WARNINGS='true'
Rscript -e 'install.packages("remotes")'
cd apis/R
Rscript -e 'remotes::install_deps(dependencies = TRUE)'
cd ../..
