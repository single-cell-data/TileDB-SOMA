#!/bin/sh
# Context: .github/workflows/*.yaml
export R_REMOTES_NO_ERRORS_FROM_WARNINGS='true'
Rscript -e 'install.packages(c("remotes", "pkgdown"))'
cd apis/R
Rscript -e 'remotes::install_deps(dependencies = TRUE)'
Rscript -e 'remotes::install_local()'
cd ../..
