# Custom vcpkg Triplets for TileDB-SOMA (CI Only)

This directory is used by CI to store dynamically generated vcpkg triplets for building release wheels.

## Purpose

Custom triplets are **not stored in the repository**. Instead, they are generated on-the-fly during CI builds in `.github/workflows/python-packaging.yml`.

## Triplets Generated (CI only)

During wheel builds, CI generates these triplets:

- `x64-linux-release.cmake` - Linux x86_64 release build
- `arm64-linux-release.cmake` - Linux ARM64 release build
- `x64-osx-release.cmake` - macOS x86_64 (Intel) release build
- `arm64-osx-release.cmake` - macOS ARM64 (Apple Silicon) release build
- `x64-windows-release.cmake` - Windows x64 release build

All triplets are configured with:

- **Build Type**: Release only (faster builds, smaller artifacts)
- **CRT Linkage**: Dynamic (runtime library linked dynamically)
- **Library Linkage**: Dynamic (shared libraries)

## Local Development

For local builds, you **don't need custom triplets**. The build will use vcpkg's default triplets automatically.

If you want to test with custom release triplets locally, you can manually create this directory and add triplet files. They will be automatically detected by the build scripts if present.

## CI Integration

The `.github/workflows/python-packaging.yml` workflow creates these triplets before running `cibuildwheel`. The `before-build` hooks in `pyproject.toml` check for this directory and copy the triplets to vcpkg if they exist.
