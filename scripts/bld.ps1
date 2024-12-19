# A script to build and install the C++ library.

<#
.SYNOPSIS
This is a Powershell script to build libtiledbsoma on Windows.

.DESCRIPTION
This script will check dependencies, and run the CMake build generator
to generate a Visual Studio solution file for libtiledbsoma, and builds
it and installs it.

.PARAMETER Prefix
Installs files in tree rooted at PREFIX (defaults to dist).

.PARAMETER Configuration
Configuration to build libtiledbsoma (defaults to Release).

.PARAMETER TileDBLocation
Location of TileDB installation (defaults to downloading it from GitHub Releases).

.LINK
https://github.com/single-cell-data/TileDB-SOMA

.EXAMPLE
bld.ps1 -Prefix "\path\to\install" -Configuration Release

#>

[CmdletBinding()]
Param(
    [string]$Configuration = 'Release',
    [string]$Prefix = '',
    [string]$TileDBLocation = '',
    [string]$TileDBRemoveDeprecations = '',
    [Alias('J')]
    [int]
    $BuildProcesses = $env:NUMBER_OF_PROCESSORS
)

$ErrorActionPreference = 'Stop'

$ExtraOpts = ''

if ($Configuration -eq 'Debug' -and $TileDBLocation -ne '') {
    $ExtraOpts += ' --debug'
}

if ($Prefix -ne '') {
    $ExtraOpts += " -DCMAKE_INSTALL_PREFIX=$Prefix -DOVERRIDE_INSTALL_PREFIX=OFF"
} else {
    $ExtraOpts += " -DOVERRIDE_INSTALL_PREFIX=ON"
}

if ($TileDBLocation -ne '') {
    $ExtraOpts += " -DTileDB_DIR=$TileDBLocation -DFORCE_BUILD_TILEDB=OFF"
}

if ($RemoveTileDBDeprecated -ne '') {
    $ExtraOpts += " -DTILEDB_REMOVE_DEPRECATIONS=$TileDBRemoveDeprecations"
}

Write-Host "Building $Configuration build"

$BuildDir = (Get-Item -Path $PSScriptRoot).Parent.FullName + '\build'
$SourceDir = (Get-Item -Path $PSScriptRoot).Parent.FullName + '\libtiledbsoma'

Remove-Item -Recurse -Force -ErrorAction SilentlyContinue $BuildDir

cmake -B $BuildDir -S $SourceDir -DCMAKE_BUILD_TYPE=$Configuration $ExtraOpts
cmake --build $BuildDir --config $Configuration -j $BuildProcesses
cmake --build $BuildDir --config $Configuration --target install-libtiledbsoma
cmake --build $BuildDir/libtiledbsoma --config $Configuration --target build_tests -j $BuildProcesses
