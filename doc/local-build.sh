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

arg() {
  echo "$1" | sed "s/^${2-[^=]*=}//" | sed "s/:/;/g"
}

# Display bootstrap usage
usage() {
echo '
Usage: '"$0"' [<options>]
Options: [defaults in brackets after descriptions]
Configuration:
    --help                          print this message
'
    exit 10
}

# Parse arguments
while test $# != 0; do
    case "$1" in
    --help) usage ;;
    *) die "Unknown option: $1" ;;
    esac
    shift
done

setup_venv() {
  if [ ! -d "${venv_dir}" ]; then
    python -m virtualenv "${venv_dir}" || die "could not create virtualenv"
  fi
  source "${venv_dir}/bin/activate" || die "could not activate virtualenv"
  pip install 'Sphinx==1.6.7' \
       'breathe' \
       'sphinx_rtd_theme' \
       -r requirements_doc.txt || die "could not install doc dependencies"
}

build_ext() {
    pushd "${ext_dir}"
    pip install apis/python || die "could not install tiledbsoma-py"
    popd
}

build_site() {
  sphinx-build -E -T -b html -d ${build_dir}/doctrees -D language=en ${source_dir} ${build_dir}/html || \
      die "could not build sphinx site"
}

run() {
  setup_venv
  build_ext
  build_site
  echo "Build complete. Open '${build_dir}/html/index.html' in your browser."
}

run
