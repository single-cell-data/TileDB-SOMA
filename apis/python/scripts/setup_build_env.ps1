# Setup script for vcpkg build environment (Windows)
# This script bootstraps vcpkg and sets up environment variables

$ErrorActionPreference = "Stop"

# Get the repository root (assuming this script is in apis/python/scripts/)
$ScriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
$PythonDir = Split-Path -Parent $ScriptDir
$RepoRoot = Split-Path -Parent (Split-Path -Parent $PythonDir)

# Check if vcpkg is already set up
if ($env:VCPKG_ROOT) {
    Write-Host "Using existing VCPKG_ROOT: $env:VCPKG_ROOT"
    $env:CMAKE_TOOLCHAIN_FILE = "$env:VCPKG_ROOT\scripts\buildsystems\vcpkg.cmake"
} elseif ($env:CMAKE_TOOLCHAIN_FILE) {
    Write-Host "Using existing CMAKE_TOOLCHAIN_FILE: $env:CMAKE_TOOLCHAIN_FILE"
    $vcpkgPath = Split-Path -Parent (Split-Path -Parent $env:CMAKE_TOOLCHAIN_FILE)
    $env:VCPKG_ROOT = $vcpkgPath
} else {
    # Try to find vcpkg in common locations
    $vcpkgRoot = $null
    if (Test-Path "$RepoRoot\vcpkg") {
        $vcpkgRoot = "$RepoRoot\vcpkg"
    } elseif (Test-Path "$env:USERPROFILE\vcpkg") {
        $vcpkgRoot = "$env:USERPROFILE\vcpkg"
    } else {
        Write-Error "vcpkg not found. Please set VCPKG_ROOT or CMAKE_TOOLCHAIN_FILE environment variable."
        exit 1
    }
    
    $env:VCPKG_ROOT = $vcpkgRoot
    
    # Bootstrap vcpkg if needed
    if (-not (Test-Path "$vcpkgRoot\vcpkg.exe")) {
        Write-Host "Bootstrapping vcpkg..."
        & "$vcpkgRoot\bootstrap-vcpkg.bat"
    }
    
    $env:CMAKE_TOOLCHAIN_FILE = "$vcpkgRoot\scripts\buildsystems\vcpkg.cmake"
}

# Set vcpkg triplet for Windows
if (-not $env:VCPKG_TARGET_TRIPLET) {
    $env:VCPKG_TARGET_TRIPLET = "x64-windows-release"
}

Write-Host "VCPKG_ROOT: $env:VCPKG_ROOT"
Write-Host "CMAKE_TOOLCHAIN_FILE: $env:CMAKE_TOOLCHAIN_FILE"
Write-Host "VCPKG_TARGET_TRIPLET: $env:VCPKG_TARGET_TRIPLET"

