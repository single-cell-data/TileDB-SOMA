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
ENV PATH="/.venv/bin:$PATH" VCPKG_FORCE_SYSTEM_BINARIES=1

ARG v=''
RUN if [ -z $v ]; then pip install tiledbsoma; else pip install tiledbsoma==$v; fi
ENTRYPOINT [ "python", "-c", "import tiledbsoma; tiledbsoma.show_package_versions()" ]
