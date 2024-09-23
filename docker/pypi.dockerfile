ARG FROM=ubuntu:24.04
FROM $FROM

ENV DEBIAN_FRONTEND=noninteractive

RUN apt update -y \
 && apt install -y software-properties-common \
 && add-apt-repository -y ppa:deadsnakes/ppa \
 && apt update -y \
 && apt install -y cmake curl g++ git make ninja-build pkg-config tar unzip zip \
 && apt install -y python3.12-dev python3.12-venv python3-pip \
 && update-alternatives --install /usr/bin/python python /usr/bin/python3.12 1 \
 && python -m venv .venv
ENV PATH="/.venv/bin:$PATH" VCPKG_FORCE_SYSTEM_BINARIES=1

ARG v=""
RUN if [ -z "$v" ]; then pip install tiledbsoma; else pip install tiledbsoma==$v; fi

ENTRYPOINT [ "python", "-c", "import tiledbsoma; tiledbsoma.show_package_versions()" ]
