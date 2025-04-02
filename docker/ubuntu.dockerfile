# Example Ubuntu 24.04 image with `tiledbsoma` installed from PyPI
ARG ubuntu=24.04
FROM ubuntu:$ubuntu

ENV DEBIAN_FRONTEND=noninteractive

ARG python=3.12
RUN apt update -y \
 && apt install -y software-properties-common \
 && add-apt-repository -y ppa:deadsnakes/ppa \
 && apt update -y \
 && apt install -y cmake curl g++ git make ninja-build pkg-config tar unzip zip \
 && apt install -y python${python}-dev python${python}-venv python3-pip \
 && update-alternatives --install /usr/bin/python python /usr/bin/python${python} 1 \
 && python -m venv .venv
ENV PATH="/.venv/bin:$PATH"

# Require CMake<4 on ARM Linux: https://github.com/single-cell-data/TileDB-SOMA/issues/3890
RUN if [ "$(uname -m)" = aarch64 ]; then echo 'cmake<4' > constraints.txt; else touch constraints.txt; fi
ENV PIP_CONSTRAINT=/constraints.txt
# vcpkg requires this on ARM, elsewhere it's a no-op (as `pip install tiledbsoma` below installs wheels, skips vcpkg bootstrap)
ENV VCPKG_FORCE_SYSTEM_BINARIES=1

ARG v=''
RUN if [ -z $v ]; then pip install tiledbsoma; else pip install tiledbsoma==$v; fi
ENTRYPOINT [ "python", "-c", "import tiledbsoma; tiledbsoma.show_package_versions()" ]
