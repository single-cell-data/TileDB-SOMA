#
# CMakeLists.txt
#
#
# The MIT License
#
# Copyright (c) 2022 TileDB, Inc.
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
#

############################################################
# CMake setup
############################################################

cmake_minimum_required(VERSION 3.21)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules/")

# TileDBSC Options
option(TILEDBSC_ENABLE_PYTHON "Enables Python API bindings" ON)

# Superbuild option must be on by default.
option(SUPERBUILD "If true, perform a superbuild (builds all missing dependencies)." ON)
option(CMAKE_IDE "(Used for CLion builds). Disables superbuild and sets the EP install dir." OFF)
option(FORCE_BUILD_TILEDB "Forces a local build of TileDB instead of searching system paths." ON)
option(TILEDB_FORCE_ALL_DEPS "Forces a local build of TileDB instead of searching system paths." OFF)
option(DOWNLOAD_TILEDB_PREBUILT "If tiledb is being super built, this controls downloading prebuilt artifacts or building from source" ON)

# TileDB BUILD Options
option(TILEDB_S3 "Enables S3/minio support using aws-cpp-sdk" ON)
option(TILEDB_AZURE "Enables Azure Storage support using azure-storage-cpp" ON)
option(TILEDB_GCS "Enables GCS Storage support using google-cloud-cpp" OFF)
option(TILEDB_HDFS "Enables HDFS support using the official Hadoop JNI bindings" OFF)
option(TILEDB_WERROR "Enables the -Werror flag during compilation." OFF)
option(TILEDB_SERIALIZATION "If true, enables building with support for query serialization" OFF)
option(OVERRIDE_INSTALL_PREFIX "Ignores the setting of CMAKE_INSTALL_PREFIX and sets a default prefix" ON)
option(ENABLE_ARROW_EXPORT "Installs an extra header for exporting in-memory results with Apache Arrow" ON)
option(TILEDB_LOG_OUTPUT_ON_FAILURE "If true, print error logs if dependency sub-project build fails" ON)
option(TILEDBSC_TESTING "Enable tests" ON)

# Set C++20 as required standard for all C++ targets (C++17 minimum is required to use the TileDB C++ API).
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF) # Don't use GNU extensions

# Build with fPIC
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Release builds by default.
if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RELEASE)
endif()

# Use @rpath on macOS for building shared libraries.
if (APPLE)
  set(CMAKE_MACOSX_RPATH ON)

  # Don't allow macOS .frameworks to be used for dependencies.
  set(CMAKE_FIND_FRAMEWORK NEVER)
endif()

# Set -fvisibility=hidden (or equivalent) flags by default.
set(CMAKE_C_VISIBILITY_PRESET hidden)
set(CMAKE_CXX_VISIBILITY_PRESET hidden)

# Root directory default installation prefix
if (OVERRIDE_INSTALL_PREFIX)
  set(PREFIX_REL_PATH "${CMAKE_SOURCE_DIR}/../dist")
  get_filename_component(DEFAULT_PREFIX "${PREFIX_REL_PATH}" ABSOLUTE)
  set(CMAKE_INSTALL_PREFIX "${DEFAULT_PREFIX}" CACHE PATH "Default install prefix" FORCE)
  message(STATUS "Using default install prefix ${CMAKE_INSTALL_PREFIX}. To control CMAKE_INSTALL_PREFIX, set OVERRIDE_INSTALL_PREFIX=OFF")
endif()
message(STATUS "Install prefix is ${CMAKE_INSTALL_PREFIX}.")

if (FORCE_BUILD_TILEDB)
  message(STATUS "Skipping search for TileDB, building it as an external project. To use system TileDB, set FORCE_BUILD_TILEDB=OFF")
endif()


############################################################
# Superbuild setup
############################################################

# Search the externals install directory for dependencies.
list(APPEND CMAKE_PREFIX_PATH "${EP_INSTALL_PREFIX}")

# If this is an in-IDE build, we need to disable the superbuild and explicitly
# set the EP base dir. The normal 'cmake && make' process won't need this step,
# it is for better CLion support of the superbuild architecture.
if (CMAKE_IDE)
  set(SUPERBUILD OFF)
  set(EP_BASE "${CMAKE_CURRENT_BINARY_DIR}/externals")
endif()

if (SUPERBUILD)
  project(TileDB-SC-Superbuild)
  message(STATUS "Starting TileDB-SingleCell superbuild.")
  include("cmake/Superbuild.cmake")
  return()
endif()

project(TileDB-SC)
message(STATUS "Starting TileDB-SC regular build.")

# Paths to locate the installed external projects.
set(EP_SOURCE_DIR "${EP_BASE}/src")
set(EP_INSTALL_PREFIX "${EP_BASE}/install")

############################################################
# Compile options/definitions for all targets
############################################################

# Set compiler flags
if (MSVC)
  # We disable some warnings that are not present in gcc/clang -Wall:
  #   C4101: unreferenced local variable
  #   C4146: unary minus operator applied to unsigned type
  #   C4244: conversion warning of floating point to integer type.
  #   C4251: C++ export warning
  #   C4456: local variable hiding previous local variable
  #   C4457: local variable hiding function parameter
  #   C4702: unreachable code
  #   C4800: warning implicit cast int to bool
  #   C4996: deprecation warning about e.g. sscanf.
  add_compile_options(/W4 /wd4101 /wd4146 /wd4244 /wd4251 /wd4456 /wd4457 /wd4702 /wd4800 /wd4996)
  # Warnings as errors:
  if (TILEDB_WERROR)
    add_compile_options(/WX)
  endif()
  # Disable GDI (which we don't need, and causes some macro
  # re-definition issues if wingdi.h is included)
  add_compile_options(/DNOGDI)
  # Add /MPn flag from CMake invocation (if defined).
  add_compile_options(${MSVC_MP_FLAG})
  # Build-specific flags
  add_compile_options(
                      "$<$<CONFIG:Debug>:/DDEBUG /Od /Zi /bigobj>"
                      "$<$<CONFIG:Release>:/DNDEBUG /Ox>"
                      "$<$<CONFIG:RelWithDebInfo>:/DNDEBUG /Ox /Zi>"
                      )
else()
  add_compile_options(-Wall -Wextra)
  if (TILEDB_WERROR)
    add_compile_options(-Werror)
  endif()
  # Build-specific flags
  if (CMAKE_BUILD_TYPE MATCHES "Debug")
    add_compile_options(-DDEBUG -O0 -g3 -ggdb3 -gdwarf-3)
  elseif (CMAKE_BUILD_TYPE MATCHES "Release")
    add_compile_options(-DNDEBUG -O3)
  elseif (CMAKE_BUILD_TYPE MATCHES "RelWithDebInfo")
    add_compile_options(-DNDEBUG -O3 -g3 -ggdb3 -gdwarf-3)
  elseif (CMAKE_BUILD_TYPE MATCHES "Coverage")
    add_compile_options(-DDEBUG -g3 -gdwarf-3 --coverage)
  endif()

  # Use -Wno-literal-suffix on Linux with C++ sources.
  if (NOT CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    add_compile_options($<$<COMPILE_LANGUAGE:CXX>:-Wno-literal-suffix>)
  endif()
endif()

# Definitions for all targets
add_definitions(-D_FILE_OFFSET_BITS=64)

# AVX2 flag
include(CheckAVX2Support)
CheckAVX2Support()
if (COMPILER_SUPPORTS_AVX2)
  add_compile_options(${COMPILER_AVX2_FLAG})
endif()

############################################################
# Regular build
############################################################

# Enable testing
enable_testing()

add_subdirectory(src)
add_subdirectory(test)
