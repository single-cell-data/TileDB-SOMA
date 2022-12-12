#!/bin/bash

set -euxo pipefail

PKG=tiledbsoma

# first build the sdist
cd "$(dirname "$0")"
rm -rf dist/"${PKG}"*
python3 setup.py sdist

# prepare manylinux docker image with build dependencies
docker pull quay.io/pypa/manylinux2014_x86_64
echo '
FROM quay.io/pypa/manylinux2014_x86_64
RUN yum install -y python3-wheel python3-pip
RUN pip3 install "pybind11[global]"
' | docker build -t "${PKG}_manylinux_builder" -

# use manylinux container to build wheel and copy it back out to dist/
docker run --rm -it -u "$(id -u):$(id -g)" -v "$(pwd)/dist":/dist "${PKG}_manylinux_builder" bash -c "
set -euxo pipefail
mkdir /tmp/bdist
cd /tmp/bdist
tar zvxf /dist/${PKG}-*.tar.gz
cd ${PKG}*
python3 setup.py bdist_wheel
cp dist/${PKG}-*.whl /dist
"
