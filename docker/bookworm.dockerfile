# Debian 12 image with GCC 13
ARG gcc=13
FROM gcc:$gcc-bookworm

ENV DEBIAN_FRONTEND=noninteractive VCPKG_FORCE_SYSTEM_BINARIES=1
RUN apt update -y \
 && apt install -y cmake curl g++ git make ninja-build pkg-config tar unzip zip \
 && apt install -y python3 python-is-python3 python3-pip python3-venv \
 && python -m venv .venv
ENV PATH=/.venv/bin:$PATH

# tiledbsoma>=1.16 requires GCC 13
ARG v=''
RUN if [ -z $v ]; then pip install tiledbsoma; else pip install tiledbsoma==$v; fi
ENTRYPOINT [ "python", "-c", "import tiledbsoma; tiledbsoma.show_package_versions()" ]
