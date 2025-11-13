#!/usr/bin/env bash
# Setup script for vcpkg build environment (Linux/macOS)
# This script bootstraps vcpkg and sets up environment variables

set -eu -o pipefail

# Get the repository root (assuming this script is in apis/python/scripts/)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PYTHON_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
REPO_ROOT="$(cd "${PYTHON_DIR}/../.." && pwd)"

# Check if vcpkg is already set up
if [ -n "${VCPKG_ROOT:-}" ]; then
    echo "Using existing VCPKG_ROOT: ${VCPKG_ROOT}"
    export CMAKE_TOOLCHAIN_FILE="${VCPKG_ROOT}/scripts/buildsystems/vcpkg.cmake"
elif [ -n "${CMAKE_TOOLCHAIN_FILE:-}" ]; then
    echo "Using existing CMAKE_TOOLCHAIN_FILE: ${CMAKE_TOOLCHAIN_FILE}"
    export VCPKG_ROOT="$(dirname "$(dirname "${CMAKE_TOOLCHAIN_FILE}")")"
else
    # Try to find vcpkg in common locations
    if [ -d "${REPO_ROOT}/vcpkg" ]; then
        export VCPKG_ROOT="${REPO_ROOT}/vcpkg"
    elif [ -d "${HOME}/vcpkg" ]; then
        export VCPKG_ROOT="${HOME}/vcpkg"
    else
        echo "Error: vcpkg not found. Please set VCPKG_ROOT or CMAKE_TOOLCHAIN_FILE environment variable."
        exit 1
    fi
    
    # Bootstrap vcpkg if needed
    if [ ! -f "${VCPKG_ROOT}/vcpkg" ] && [ ! -f "${VCPKG_ROOT}/vcpkg.exe" ]; then
        echo "Bootstrapping vcpkg..."
        if [ "$(uname)" == "Darwin" ]; then
            "${VCPKG_ROOT}/bootstrap-vcpkg.sh"
        else
            "${VCPKG_ROOT}/bootstrap-vcpkg.sh"
        fi
    fi
    
    export CMAKE_TOOLCHAIN_FILE="${VCPKG_ROOT}/scripts/buildsystems/vcpkg.cmake"
fi

# Set vcpkg triplet based on platform and architecture
if [ -z "${VCPKG_TARGET_TRIPLET:-}" ]; then
    if [ "$(uname)" == "Darwin" ]; then
        arch=$(uname -m)
        if [ "$arch" == "arm64" ]; then
            export VCPKG_TARGET_TRIPLET="arm64-osx-release"
        else
            export VCPKG_TARGET_TRIPLET="x64-osx-release"
        fi
    elif [ "$(uname)" == "Linux" ]; then
        arch=$(uname -m)
        if [ "$arch" == "aarch64" ]; then
            export VCPKG_TARGET_TRIPLET="arm64-linux-release"
        else
            export VCPKG_TARGET_TRIPLET="x64-linux-release"
        fi
    else
        echo "Warning: Unknown platform, using default triplet"
        export VCPKG_TARGET_TRIPLET="x64-linux-release"
    fi
fi

# Set macOS architecture if on macOS
if [ "$(uname)" == "Darwin" ] && [ -z "${CMAKE_OSX_ARCHITECTURES:-}" ]; then
    export CMAKE_OSX_ARCHITECTURES="$(uname -m)"
fi

echo "VCPKG_ROOT: ${VCPKG_ROOT}"
echo "CMAKE_TOOLCHAIN_FILE: ${CMAKE_TOOLCHAIN_FILE}"
echo "VCPKG_TARGET_TRIPLET: ${VCPKG_TARGET_TRIPLET}"
if [ "$(uname)" == "Darwin" ]; then
    echo "CMAKE_OSX_ARCHITECTURES: ${CMAKE_OSX_ARCHITECTURES}"
fi

