# Setup build environment for cibuildwheel on Windows
# This script bootstraps vcpkg and sets up the CMAKE_TOOLCHAIN_FILE

Write-Host "=== Setting up build environment ===" -ForegroundColor Green

# Determine the repository root (2 levels up from this script)
$ScriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
$RepoRoot = Split-Path -Parent (Split-Path -Parent $ScriptDir)
$VcpkgRoot = Join-Path $RepoRoot "vcpkg"

Write-Host "Script directory: $ScriptDir"
Write-Host "Repository root: $RepoRoot"
Write-Host "vcpkg root: $VcpkgRoot"

# Check if vcpkg exists
if (-not (Test-Path $VcpkgRoot)) {
    Write-Host "Cloning vcpkg..." -ForegroundColor Yellow
    git clone https://github.com/microsoft/vcpkg.git $VcpkgRoot
}

# Bootstrap vcpkg if not already done
$VcpkgExe = Join-Path $VcpkgRoot "vcpkg.exe"
if (-not (Test-Path $VcpkgExe)) {
    Write-Host "Bootstrapping vcpkg..." -ForegroundColor Yellow
    $BootstrapScript = Join-Path $VcpkgRoot "bootstrap-vcpkg.bat"
    & $BootstrapScript
}
else {
    Write-Host "vcpkg already bootstrapped" -ForegroundColor Green
}

# Set the CMAKE_TOOLCHAIN_FILE environment variable
$ToolchainFile = Join-Path $VcpkgRoot "scripts\buildsystems\vcpkg.cmake"
if (Test-Path $ToolchainFile) {
    Write-Host "Setting CMAKE_TOOLCHAIN_FILE=$ToolchainFile" -ForegroundColor Green
    $env:CMAKE_TOOLCHAIN_FILE = $ToolchainFile
    [System.Environment]::SetEnvironmentVariable("CMAKE_TOOLCHAIN_FILE", $ToolchainFile, [System.EnvironmentVariableTarget]::Process)
    
    # Write to GitHub environment if in CI
    if ($env:GITHUB_ENV) {
        "CMAKE_TOOLCHAIN_FILE=$ToolchainFile" | Out-File -FilePath $env:GITHUB_ENV -Encoding utf8 -Append
    }
}
else {
    Write-Host "Warning: vcpkg toolchain file not found at $ToolchainFile" -ForegroundColor Red
}

Write-Host "=== Build environment setup complete ===" -ForegroundColor Green

