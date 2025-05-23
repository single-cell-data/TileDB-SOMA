# Example Debian 12 image, comes with GCC 12 and Python 3.12
FROM langchain/langgraph-api:3.12

ENV DEBIAN_FRONTEND=noninteractive
RUN apt update -y \
 && apt install -y cmake curl g++ git make ninja-build pkg-config tar unzip zip

# Require CMake<4 on ARM Linux: https://github.com/single-cell-data/TileDB-SOMA/issues/3890
RUN if [ "$(uname -m)" = aarch64 ]; then echo 'cmake<4' > constraints.txt; else touch constraints.txt; fi
ENV PIP_CONSTRAINT=/api/constraints.txt
# vcpkg requires this on ARM, elsewhere it's a no-op (as `pip install tiledbsoma` below installs wheels, skips vcpkg bootstrap)
ENV VCPKG_FORCE_SYSTEM_BINARIES=1

# tiledbsoma>=1.16 requires GCC 13, which is more involved to install on this base image; see bookworm-pypi for an
# example that works with both GCC 12 and 13
ARG version=1.15.7
RUN pip install tiledbsoma==$version
ENTRYPOINT [ "python", "-c", "import tiledbsoma; tiledbsoma.show_package_versions()" ]
