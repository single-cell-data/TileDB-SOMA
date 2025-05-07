# Example Ubuntu 24.04 image with `tiledbsoma` installed from source
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
ENV PATH="/.venv/bin:$PATH" VCPKG_FORCE_SYSTEM_BINARIES=1

ARG ref=main
RUN git clone -b "$ref" https://github.com/single-cell-data/TileDB-SOMA
WORKDIR TileDB-SOMA

# Release or Debug
ARG build=Release
RUN if ! make install build=$build; then cat /TileDB-SOMA/build/externals/src/ep_tiledb-stamp/ep_tiledb-configure-err.log; exit 1; fi

ENTRYPOINT [ "python", "scripts/show-versions.py" ]
