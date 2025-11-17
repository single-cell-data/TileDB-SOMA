#!/usr/bin/env bash
# Setup script for vcpkg build environment (Linux/macOS)
# This script bootstraps vcpkg and sets up environment variables
# Usage: source this script or eval "$(./setup_build_env.sh -x)"

set -eu -o pipefail

# Get the repository root (assuming this script is in apis/python/scripts/)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PYTHON_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
REPO_ROOT="$(cd "${PYTHON_DIR}/../.." && pwd)"

# Check if vcpkg is already set up
if [ -n "${CMAKE_TOOLCHAIN_FILE:-}" ] && [ -f "${CMAKE_TOOLCHAIN_FILE}" ]; then
    # CMAKE_TOOLCHAIN_FILE is already set and exists - use it
    echo "Using existing CMAKE_TOOLCHAIN_FILE: ${CMAKE_TOOLCHAIN_FILE}"
    VCPKG_ROOT="$(dirname "$(dirname "${CMAKE_TOOLCHAIN_FILE}")")"
elif [ -n "${VCPKG_ROOT:-}" ] && [ -d "${VCPKG_ROOT}" ]; then
    # VCPKG_ROOT is set and exists - use it
    echo "Using existing VCPKG_ROOT: ${VCPKG_ROOT}"
    CMAKE_TOOLCHAIN_FILE="${VCPKG_ROOT}/scripts/buildsystems/vcpkg.cmake"
else
    # Try to find vcpkg in common locations
    if [ -d "${REPO_ROOT}/vcpkg" ]; then
        VCPKG_ROOT="${REPO_ROOT}/vcpkg"
    elif [ -d "${PYTHON_DIR}/../vcpkg" ]; then
        # For cibuildwheel: vcpkg is cloned one level up from apis/python
        VCPKG_ROOT="${PYTHON_DIR}/../vcpkg"
    elif [ -d "${HOME}/vcpkg" ]; then
        VCPKG_ROOT="${HOME}/vcpkg"
    else
        echo "Error: vcpkg not found. Please set VCPKG_ROOT or CMAKE_TOOLCHAIN_FILE environment variable."
        echo "Searched locations:"
        echo "  - ${REPO_ROOT}/vcpkg"
        echo "  - ${PYTHON_DIR}/../vcpkg"
        echo "  - ${HOME}/vcpkg"
        echo ""
        echo "To clone vcpkg, run:"
        echo "  git clone https://github.com/microsoft/vcpkg.git ${REPO_ROOT}/vcpkg"
        exit 1
    fi
    
    # Bootstrap vcpkg if needed
    if [ ! -f "${VCPKG_ROOT}/vcpkg" ] && [ ! -f "${VCPKG_ROOT}/vcpkg.exe" ]; then
        echo "Bootstrapping vcpkg..."
        "${VCPKG_ROOT}/bootstrap-vcpkg.sh" || {
            echo "Error: Failed to bootstrap vcpkg"
            exit 1
        }
    fi
    
    # Copy custom triplets if they exist in the Python package
    if [ -d "${PYTHON_DIR}/vcpkg_triplets" ] && [ -d "${VCPKG_ROOT}/triplets" ]; then
        echo "Installing custom vcpkg triplets..."
        cp -f "${PYTHON_DIR}/vcpkg_triplets/"*-release.cmake "${VCPKG_ROOT}/triplets/" 2>/dev/null || true
    fi
    
    CMAKE_TOOLCHAIN_FILE="${VCPKG_ROOT}/scripts/buildsystems/vcpkg.cmake"
fi

# Set vcpkg triplet based on platform and architecture
# For local builds, we use standard vcpkg triplets (not the custom release ones used in CI)
# This can be overridden by setting VCPKG_TARGET_TRIPLET environment variable
if [ -z "${VCPKG_TARGET_TRIPLET:-}" ]; then
    if [ "$(uname)" == "Darwin" ]; then
        arch=$(uname -m)
        if [ "$arch" == "arm64" ]; then
            VCPKG_TARGET_TRIPLET="arm64-osx"
        else
            VCPKG_TARGET_TRIPLET="x64-osx"
        fi
    elif [ "$(uname)" == "Linux" ]; then
        arch=$(uname -m)
        if [ "$arch" == "aarch64" ]; then
            VCPKG_TARGET_TRIPLET="arm64-linux"
        else
            VCPKG_TARGET_TRIPLET="x64-linux"
        fi
    else
        echo "Warning: Unknown platform, using default triplet"
        VCPKG_TARGET_TRIPLET="x64-linux"
    fi
fi

# Set macOS architecture if on macOS
if [ "$(uname)" == "Darwin" ] && [ -z "${CMAKE_OSX_ARCHITECTURES:-}" ]; then
    CMAKE_OSX_ARCHITECTURES="$(uname -m)"
fi

# Export variables
export VCPKG_ROOT
export CMAKE_TOOLCHAIN_FILE
export VCPKG_TARGET_TRIPLET
if [ "$(uname)" == "Darwin" ]; then
    export CMAKE_OSX_ARCHITECTURES
fi

# If -x flag is passed, output export commands for eval
# This allows the script to be used with: eval "$(./setup_build_env.sh -x)"
if [ "${1:-}" = "-x" ]; then
    echo "export VCPKG_ROOT=\"${VCPKG_ROOT}\""
    echo "export CMAKE_TOOLCHAIN_FILE=\"${CMAKE_TOOLCHAIN_FILE}\""
    echo "export VCPKG_TARGET_TRIPLET=\"${VCPKG_TARGET_TRIPLET}\""
    if [ "$(uname)" == "Darwin" ] && [ -n "${CMAKE_OSX_ARCHITECTURES:-}" ]; then
        echo "export CMAKE_OSX_ARCHITECTURES=\"${CMAKE_OSX_ARCHITECTURES}\""
    fi
    exit 0
fi

echo "VCPKG_ROOT: ${VCPKG_ROOT}"
echo "CMAKE_TOOLCHAIN_FILE: ${CMAKE_TOOLCHAIN_FILE}"
echo "VCPKG_TARGET_TRIPLET: ${VCPKG_TARGET_TRIPLET}"
if [ "$(uname)" == "Darwin" ]; then
    echo "CMAKE_OSX_ARCHITECTURES: ${CMAKE_OSX_ARCHITECTURES}"
fi

