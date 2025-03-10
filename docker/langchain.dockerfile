# Example Debian 12 image, comes with Python 3.12
FROM langchain/langgraph-api:3.12
ENV DEBIAN_FRONTEND=noninteractive VCPKG_FORCE_SYSTEM_BINARIES=1
RUN apt update -y && apt install -y cmake curl g++ git make ninja-build pkg-config tar unzip zip
RUN pip install tiledbsoma
ENTRYPOINT [ "python", "-c", "import tiledbsoma; tiledbsoma.show_package_versions()" ]
