#!/bin/bash
# Setup build environment for cibuildwheel
# This script bootstraps vcpkg and sets up the CMAKE_TOOLCHAIN_FILE

set -e

echo "=== Setting up build environment ==="

# Determine the repository root (2 levels up from this script)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
VCPKG_ROOT="$REPO_ROOT/vcpkg"

echo "Script directory: $SCRIPT_DIR"
echo "Repository root: $REPO_ROOT"
echo "vcpkg root: $VCPKG_ROOT"

# Check if vcpkg exists
if [ ! -d "$VCPKG_ROOT" ]; then
    echo "Cloning vcpkg..."
    git clone https://github.com/microsoft/vcpkg.git "$VCPKG_ROOT"
fi

# Bootstrap vcpkg if not already done
if [ ! -f "$VCPKG_ROOT/vcpkg" ]; then
    echo "Bootstrapping vcpkg..."
    cd "$VCPKG_ROOT"
    ./bootstrap-vcpkg.sh
    cd -
else
    echo "vcpkg already bootstrapped"
fi

# Set the CMAKE_TOOLCHAIN_FILE environment variable
TOOLCHAIN_FILE="$VCPKG_ROOT/scripts/buildsystems/vcpkg.cmake"
if [ -f "$TOOLCHAIN_FILE" ]; then
    echo "Setting CMAKE_TOOLCHAIN_FILE=$TOOLCHAIN_FILE"
    export CMAKE_TOOLCHAIN_FILE="$TOOLCHAIN_FILE"
    
    # Write to GitHub environment if in CI
    if [ -n "$GITHUB_ENV" ]; then
        echo "CMAKE_TOOLCHAIN_FILE=$TOOLCHAIN_FILE" >> "$GITHUB_ENV"
    fi
else
    echo "Warning: vcpkg toolchain file not found at $TOOLCHAIN_FILE"
fi

echo "=== Build environment setup complete ==="

