# Debian 12 image with GCC 13
ARG gcc=13
FROM gcc:$gcc-bookworm

ENV DEBIAN_FRONTEND=noninteractive
RUN apt update -y \
 && apt install -y cmake curl g++ git make ninja-build pkg-config tar unzip zip \
 && apt install -y python3 python-is-python3 python3-pip python3-venv \
 && python -m venv .venv
ENV PATH=/.venv/bin:$PATH

# Require CMake<4 on ARM Linux: https://github.com/single-cell-data/TileDB-SOMA/issues/3890
RUN if [ "$(uname -m)" = aarch64 ]; then echo 'cmake<4' > constraints.txt; else touch constraints.txt; fi
ENV PIP_CONSTRAINT=/constraints.txt
# vcpkg requires this on ARM, elsewhere it's a no-op (as `pip install tiledbsoma` below installs wheels, skips vcpkg bootstrap)
ENV VCPKG_FORCE_SYSTEM_BINARIES=1

# tiledbsoma>=1.16 requires GCC 13
ARG v=''
RUN if [ -z $v ]; then pip install tiledbsoma; else pip install tiledbsoma==$v; fi
ENTRYPOINT [ "python", "-c", "import tiledbsoma; tiledbsoma.show_package_versions()" ]
