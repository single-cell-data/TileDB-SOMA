# Setup script for vcpkg build environment (Windows PowerShell)
# This script bootstraps vcpkg and sets up environment variables

$ErrorActionPreference = "Stop"

# Get the repository root (assuming this script is in apis/python/scripts/)
$ScriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
$PythonDir = Split-Path -Parent $ScriptDir
$RepoRoot = Split-Path -Parent (Split-Path -Parent $PythonDir)

# Check if vcpkg is already set up
if ($env:CMAKE_TOOLCHAIN_FILE -and (Test-Path $env:CMAKE_TOOLCHAIN_FILE)) {
    # CMAKE_TOOLCHAIN_FILE is already set and exists - use it
    Write-Host "Using existing CMAKE_TOOLCHAIN_FILE: $env:CMAKE_TOOLCHAIN_FILE"
    $VcpkgRoot = Split-Path -Parent (Split-Path -Parent $env:CMAKE_TOOLCHAIN_FILE)
}
elseif ($env:VCPKG_ROOT -and (Test-Path $env:VCPKG_ROOT)) {
    # VCPKG_ROOT is set and exists - use it
    Write-Host "Using existing VCPKG_ROOT: $env:VCPKG_ROOT"
    $VcpkgRoot = $env:VCPKG_ROOT
    $env:CMAKE_TOOLCHAIN_FILE = Join-Path $VcpkgRoot "scripts\buildsystems\vcpkg.cmake"
}
else {
    # Try to find vcpkg in common locations
    $VcpkgRoot = $null
    
    if (Test-Path (Join-Path $RepoRoot "vcpkg")) {
        $VcpkgRoot = Join-Path $RepoRoot "vcpkg"
    }
    elseif (Test-Path (Join-Path (Split-Path -Parent $PythonDir) "vcpkg")) {
        # For cibuildwheel: vcpkg is cloned one level up from apis/python
        $VcpkgRoot = Join-Path (Split-Path -Parent $PythonDir) "vcpkg"
    }
    elseif (Test-Path (Join-Path $env:USERPROFILE "vcpkg")) {
        $VcpkgRoot = Join-Path $env:USERPROFILE "vcpkg"
    }
    else {
        Write-Error @"
Error: vcpkg not found. Please set VCPKG_ROOT or CMAKE_TOOLCHAIN_FILE environment variable.
Searched locations:
  - $RepoRoot\vcpkg
  - $(Split-Path -Parent $PythonDir)\vcpkg
  - $env:USERPROFILE\vcpkg

To clone vcpkg, run:
  git clone https://github.com/microsoft/vcpkg.git $RepoRoot\vcpkg
"@
        exit 1
    }
    
    # Bootstrap vcpkg if needed
    $VcpkgExe = Join-Path $VcpkgRoot "vcpkg.exe"
    if (-not (Test-Path $VcpkgExe)) {
        Write-Host "Bootstrapping vcpkg..."
        $BootstrapScript = Join-Path $VcpkgRoot "bootstrap-vcpkg.bat"
        & $BootstrapScript
        if ($LASTEXITCODE -ne 0) {
            Write-Error "Error: Failed to bootstrap vcpkg"
            exit 1
        }
    }
    
    # Copy custom triplets if they exist in the Python package
    $TripletsSource = Join-Path $PythonDir "vcpkg_triplets"
    $TripletsTarget = Join-Path $VcpkgRoot "triplets"
    if ((Test-Path $TripletsSource) -and (Test-Path $TripletsTarget)) {
        Write-Host "Installing custom vcpkg triplets..."
        Copy-Item "$TripletsSource\*-release.cmake" $TripletsTarget -Force -ErrorAction SilentlyContinue
    }
    
    $env:CMAKE_TOOLCHAIN_FILE = Join-Path $VcpkgRoot "scripts\buildsystems\vcpkg.cmake"
}

# Set vcpkg triplet
# For local builds, we use standard vcpkg triplets (not the custom release ones used in CI)
if (-not $env:VCPKG_TARGET_TRIPLET) {
    # Windows is always x64 for now (use standard triplet)
    $env:VCPKG_TARGET_TRIPLET = "x64-windows"
}

# Export variables
$env:VCPKG_ROOT = $VcpkgRoot

Write-Host "VCPKG_ROOT: $env:VCPKG_ROOT"
Write-Host "CMAKE_TOOLCHAIN_FILE: $env:CMAKE_TOOLCHAIN_FILE"
Write-Host "VCPKG_TARGET_TRIPLET: $env:VCPKG_TARGET_TRIPLET"
