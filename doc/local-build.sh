#!/bin/bash

#
# Builds the ReadTheDocs documentation locally.
# Usage. Execute in this directory:
#   $ ./local-build.sh
# This creates a Python virtual env 'venv' in the current directory.
#

# Choose the default directories
source_dir="source"
build_dir="source/_build"
venv_dir="venv"
ext_dir="../apis/python"

die() {
  echo "$@" 1>&2 ; popd 2>/dev/null; exit 1
}

# Display bootstrap usage
usage() {
    echo '
    Usage: '"$0"' [-r/--reinstall] [--help]
    Options: [defaults in brackets after descriptions]
    Configuration:
        -r, --reinstall                 reinstall dependencies
        --help                          print this message
    '
    exit 10
}

# Parse arguments
reinstall=
while test $# != 0; do
    case "$1" in
    -r|--reinstall) reinstall=1 ;;
    --help) usage ;;
    *) die "Unknown option: $1" ;;
    esac
    shift
done

made_venv=
if [ ! -d "${venv_dir}" ]; then
  python -m virtualenv "${venv_dir}" || die "could not create virtualenv"
  made_venv=1
fi

source "${venv_dir}/bin/activate" || die "could not activate virtualenv"

if [ -n "${reinstall}" ] || [ -n "${made_venv}" ]; then
  pip install -r requirements_doc.txt || die "could not install doc dependencies"
  pushd "${ext_dir}"
  pip install -e . || die "could not install tiledbsoma-py"
  popd
fi

sphinx-build -E -T -b html -d ${build_dir}/doctrees -D language=en ${source_dir} ${build_dir}/html || \
  die "could not build sphinx site"

echo "Build complete. Open '${build_dir}/html/index.html' in your browser."
