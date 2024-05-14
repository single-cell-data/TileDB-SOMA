#!/bin/bash

# Builds the ReadTheDocs documentation.
# By default, creates a Python virtual env in 'venv' (if not already present).

pwd="$(pwd)"
cd "$(dirname "${BASH_SOURCE[0]}")" || exit 1
parent="$(pwd)"
root="$(dirname "$parent")"
cd "$root"

# Choose the default directories
source_dir="doc/source"
build_dir="_readthedocs"
venv_dir="doc/venv"
ext_dir="apis/python"

die() {
  echo "$@" 1>&2 ; popd 2>/dev/null; exit 1
}

usage() {
    echo '
    Usage: '"$0"' [...options | --help]
    Options:
        -b, --build-dir=DIR             build directory ['"$build_dir"']
        -o, --open                      `open` the built documentation after building
        -r, --reinstall                 reinstall dependencies
        -v, --venv-dir=DIR              virtualenv directory ['"$venv_dir"']
        -V, --no-venv                   do not create a virtualenv
        --help                          print this message
    '
    exit 10
}

# Parse arguments
open=
reinstall=
use_venv=1
while test $# != 0; do
    case "$1" in
    -b|--build-dir) shift; build_dir="$pwd/$1" ;;
    -o|--open) open=1 ;;
    -r|--reinstall) reinstall=1 ;;
    -v|--venv-dir) shift; venv_dir="$pwd/$1" ;;
    -V|--no-venv) use_venv= ;;
    --help) usage ;;
    *) die "Unknown option: $1" ;;
    esac
    shift
done

made_venv=
if [ -n "$use_venv" ]; then
    if [ -d "$venv_dir" ]; then
        echo "Using existing virtualenv in '$venv_dir'" >&2
    else
        echo "Creating virtualenv in '$venv_dir'" >&2
        python -m virtualenv "$venv_dir" || die "could not create virtualenv"
        made_venv=1
    fi
    source "$venv_dir/bin/activate" || die "could not activate virtualenv"
fi

if [ -n "$reinstall" ] || [ -n "$made_venv" ]; then
  pip install -r doc/requirements_doc.txt || die "could not install doc dependencies"
  pushd "$ext_dir"
  pip install -e . || die "could not install tiledbsoma-py"
  popd
fi

sphinx-build -E -T -b html -d "$build_dir/doctrees" -D language=en "$source_dir" "$build_dir/html" || \
  die "could not build sphinx site"

index_html="$build_dir/html/index.html"
if [ -n "$open" ]; then
  if which xdg-open > /dev/null; then
    xdg-open "$index_html"
  elif which open > /dev/null; then
    open "$index_html"
  else
    echo "Could not open the built documentation. Open '$index_html' in your browser."
  fi
else
    echo "Build complete. Open '$index_html' in your browser."
fi
