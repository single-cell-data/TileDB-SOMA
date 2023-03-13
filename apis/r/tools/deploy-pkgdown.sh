#!/bin/bash

set -euo pipefail

cd apis/r

git config --local user.name "$GITHUB_ACTOR"
git config --local user.email "$GITHUB_ACTOR@users.noreply.github.com"
Rscript -e 'pkgdown::deploy_to_branch(new_process = FALSE)'

cd ../..
