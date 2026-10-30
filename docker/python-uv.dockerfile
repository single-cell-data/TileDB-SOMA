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
# Explicit scanpy>=1.11 pin is required, otherwise `uv` gets confused and fails trying to install Scanpy 1.9 → Numba 0.53 → llvmlite 0.36 → Python 3.8
RUN uv pip install --system tiledbsoma 'scanpy>=1.11'

ENTRYPOINT [ "python", "-c", "import tiledbsoma; tiledbsoma.show_package_versions()" ]
