ARG python=3.13
FROM python:$python

# Required by vcpkg source build on ARM
ENV VCPKG_FORCE_SYSTEM_BINARIES=1
# Required by vcpkg source build
RUN apt update -y && apt install -y curl ninja-build tar unzip zip

# `tiledbsoma>=1.16` requires GCC ≥13
RUN echo "deb http://deb.debian.org/debian testing main" > /etc/apt/sources.list.d/testing.list && \
    echo "deb http://deb.debian.org/debian unstable main" > /etc/apt/sources.list.d/unstable.list && \
    echo 'APT::Default-Release "stable";' > /etc/apt/apt.conf.d/99defaultrelease && \
    apt update -y && \
    apt install -y -t unstable gcc-13 g++-13 && \
    update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-13 100 && \
    update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-13 100 && \
    update-alternatives --set gcc /usr/bin/gcc-13 && \
    update-alternatives --set g++ /usr/bin/g++-13

RUN pip install --upgrade pip \
  && pip install uv \
  && uv init example

WORKDIR example
# Required on ARM, until ARM Linux wheels are published:
# - https://github.com/single-cell-data/TileDB-SOMA/issues/3890
# - https://github.com/single-cell-data/TileDB-SOMA/issues/3909
RUN echo "cmake<4" > constraints.txt
RUN uv pip install --system --build-constraints constraints.txt tiledbsoma

ENTRYPOINT [ "python", "-c", "import tiledbsoma; tiledbsoma.show_package_versions()" ]
