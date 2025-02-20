#
# CMakeLists.txt
#
#
# The MIT License
#
# Copyright (c) 2022-2024 TileDB, Inc.
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

# ###########################################################
# CMake setup
# ###########################################################

cmake_minimum_required(VERSION 3.21)

if(CMAKE_VERSION VERSION_GREATER_EQUAL "3.24.0")
  cmake_policy(SET CMP0135 NEW)
endif()

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules/")

set(TILEDBSOMA_CMAKE_INPUTS_DIR "${CMAKE_CURRENT_SOURCE_DIR}/cmake/inputs")

# Check for and set environment variables.
if (NOT DEFINED TILEDBSOMA_ENABLE_WERROR AND DEFINED $ENV{TILEDBSSOMA_ENABLE_WERROR})
  set(TILEDBSOMA_ENABLE_WERROR "$ENV{TILEDBSOMA_ENABLE_WERROR}")
endif()

# TileDB-SOMA Options
option(TILEDBSOMA_BUILD_PYTHON "Build Python API bindings" ON)
option(TILEDBSOMA_BUILD_CLI "Build tiledbsoma CLI tool" ON)
option(TILEDBSOMA_BUILD_STATIC "Build a static library; otherwise, shared library" OFF)
option(TILEDBSOMA_ENABLE_TESTING "Enable tests" ON)
option(TILEDBSOMA_ENABLE_WERROR "Enables the -Werror flag during compilation." OFF)

# Superbuild option must be on by default.
option(SUPERBUILD "If true, perform a superbuild (builds all missing dependencies)." ON)
option(CMAKE_IDE "(Used for CLion builds). Disables superbuild and sets the EP install dir." OFF)
option(FORCE_BUILD_TILEDB "Forces a local build of TileDB instead of searching system paths." OFF)
option(TILEDB_FORCE_ALL_DEPS "Forces a local build of TileDB instead of searching system paths." OFF)
option(DOWNLOAD_TILEDB_PREBUILT "If tiledb is being super built, this controls downloading prebuilt artifacts or building from source" ON)
option(SPDLOG_LINK_SHARED "Link installed spdlog as an shared library instead of using it as header only" OFF)

# TileDB BUILD Options
# NOTE: TileDB Embedded version is controlled in cmake/Modules/FindTileDB_EP.cmake
option(TILEDB_S3 "Enables S3/minio support using aws-cpp-sdk" ON)
option(TILEDB_AZURE "Enables Azure Storage support using azure-storage-cpp" ON)
option(TILEDB_GCS "Enables GCS Storage support using google-cloud-cpp" OFF)
option(TILEDB_HDFS "Enables HDFS support using the official Hadoop JNI bindings" OFF)
option(TILEDB_WERROR "Enables the -Werror flag during compilation." OFF)
option(TILEDB_REMOVE_DEPRECATIONS "If true, do not build deprecated APIs." OFF)
option(TILEDB_SERIALIZATION "If true, enables building with support for query serialization" ON)
option(TILEDB_VERBOSE "If true, sets default logging to verbose for TileDB" OFF)
option(OVERRIDE_INSTALL_PREFIX "Ignores the setting of CMAKE_INSTALL_PREFIX and sets a default prefix" OFF)
option(ENABLE_ARROW_EXPORT "Installs an extra header for exporting in-memory results with Apache Arrow" ON)
option(TILEDB_LOG_OUTPUT_ON_FAILURE "If true, print error logs if dependency sub-project build fails" ON)
option(TILEDB_SANITIZER "Sanitizer to use in TILEDB. ")

# Enable compiler cache to speed up recompilation
find_program(CCACHE_FOUND ccache)

if(CCACHE_FOUND)
  set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE ccache)
  set_property(GLOBAL PROPERTY RULE_LAUNCH_LINK ccache)
endif()

# Set C++20 as required standard for all C++ targets (C++17 minimum is required to use the TileDB C++ API).
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF) # Don't use GNU extensions

# Build with fPIC
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Set default builds/configuration to be Release.
get_property(is_multi_config GLOBAL PROPERTY GENERATOR_IS_MULTI_CONFIG)
if (is_multi_config)
  set(CMAKE_CONFIGURATION_TYPES
      "Release;Debug;RelWithDebInfo;ASAN;TSAN;LSAN;UBSAN;MSAN"
      CACHE
      STRING
      "Semi-colon separate list of build types for multi-configuration generators."
  )
  set(CMAKE_MAP_IMPORTED_CONFIG_ASAN Debug)
  set(CMAKE_MAP_IMPORTED_CONFIG_TSAN Debug)
  set(CMAKE_MAP_IMPORTED_CONFIG_LSAN Debug)
  set(CMAKE_MAP_IMPORTED_CONFIG_UBSAN Debug)
  set(CMAKE_MAP_IMPORTED_CONFIG_MSAN Debug)
else()
    set(CMAKE_BUILD_TYPE
      "Release" CACHE STRING "Build type for single-configuration generators."
    )
endif()

# Use @rpath on macOS for building shared libraries.
if(APPLE)
  set(CMAKE_MACOSX_RPATH ON)

  # Set minimum macOS version to enable certain C++20 features
  # set(CMAKE_OSX_DEPLOYMENT_TARGET 13.3)

  # Don't allow macOS .frameworks to be used for dependencies.
  set(CMAKE_FIND_FRAMEWORK NEVER)
endif()

# Set -fvisibility=hidden (or equivalent) flags by default.
#set(CMAKE_C_VISIBILITY_PRESET hidden)
#set(CMAKE_CXX_VISIBILITY_PRESET hidden)

# Root directory default installation prefix
if(OVERRIDE_INSTALL_PREFIX OR NOT CMAKE_INSTALL_PREFIX)
  set(PREFIX_REL_PATH "${CMAKE_SOURCE_DIR}/../dist")
  get_filename_component(DEFAULT_PREFIX "${PREFIX_REL_PATH}" ABSOLUTE)
  set(CMAKE_INSTALL_PREFIX "${DEFAULT_PREFIX}" CACHE PATH "Default install prefix" FORCE)
  message(STATUS "Using default install prefix ${CMAKE_INSTALL_PREFIX}. To control CMAKE_INSTALL_PREFIX, set OVERRIDE_INSTALL_PREFIX=OFF")
endif()

message(STATUS "Install prefix is ${CMAKE_INSTALL_PREFIX}.")

if(FORCE_BUILD_TILEDB)
  message(STATUS "Skipping search for TileDB, building it as an external project. To use system TileDB, set FORCE_BUILD_TILEDB=OFF")
endif()

# Export symbols without decorating the code
set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)

# Append path to linked library files to the rpath
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH ON)

# ###########################################################
# Superbuild setup
# ###########################################################

# Search the externals install directory for dependencies.
list(APPEND CMAKE_PREFIX_PATH "${EP_INSTALL_PREFIX}")

# If this is an in-IDE build, we need to disable the superbuild and explicitly
# set the EP base dir. The normal 'cmake && make' process won't need this step,
# it is for better CLion support of the superbuild architecture.
if(CMAKE_IDE)
  set(SUPERBUILD OFF)
  set(EP_BASE "${CMAKE_CURRENT_BINARY_DIR}/externals")
endif()

if(SUPERBUILD)
  project(TileDB-SOMA-Superbuild)
  message(STATUS "Starting TileDB-SOMA superbuild.")
  include("cmake/Superbuild.cmake")
  return()
endif()

project(TileDB-SOMA)
message(STATUS "Starting TileDB-SOMA regular build.")

# Paths to locate the installed external projects.
set(EP_SOURCE_DIR "${EP_BASE}/src")
set(EP_INSTALL_PREFIX "${EP_BASE}/install")


# ###########################################################
# Check compiler version supports full C++20 standard.
# ###########################################################

set(GCC_MIN 13.0)
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS GCC_MIN)
        message(FATAL_ERROR "GNU GCC must be at least version ${GCC_MIN}. Found version ${CMAKE_CXX_COMPILER_VERSION}.")
    endif()
endif()


# ###########################################################
# Compile options/definitions for all targets
# ###########################################################

# Set compiler flags
if(MSVC)
  # We disable some warnings that are not present in gcc/clang -Wall:
  # C4101: unreferenced local variable
  # C4146: unary minus operator applied to unsigned type
  # C4244: conversion warning of floating point to integer type.
  # C4251: C++ export warning
  # C4456: local variable hiding previous local variable
  # C4457: local variable hiding function parameter
  # C4702: unreachable code
  # C4800: warning implicit cast int to bool
  # C4996: deprecation warning about e.g. sscanf.
  set(TILEDBSOMA_COMPILE_OPTIONS /W4 /wd4101 /wd4146 /wd4244 /wd4251 /wd4456 /wd4457 /wd4702 /wd4800 /wd4996)

  # Warnings as errors:
  if(${TILEDBSOMA_ENABLE_WERROR})
    set(TILEDBSOMA_WERROR_OPTION /WX)
  else()
    set(TILEDBSOMA_WERROR_OPTION "")
  endif()

  # Disable GDI (which we don't need, and causes some macro
  # re-definition issues if wingdi.h is included)
  list(APPEND TILEDBSOMA_COMPILE_OPTIONS /DNOGDI)

  # Add /MPn flag from CMake invocation (if defined).
  list(APPEND TILEDBSOMA_COMPILE_OPTIONS ${MSVC_MP_FLAG})

  # Build-specific flags: Must use generator expressions to support multi-configuration
  # build generators that set the config type at build time.
  list(APPEND
    TILEDBSOMA_COMPILE_OPTIONS
    $<$<CONFIG:Debug>:/DDEBUG /Od /Zi /bigobj>
  )
  list(APPEND
    TILEDBSOMA_COMPILE_OPTIONS
    $<$<CONFIG:Release>:/DNDEBUG /Ox>
  )
  list(APPEND
    TILEDBSOMA_COMPILE_OPTIONS
    $<$<CONFIG:RelWithDebInfo>:/DNDEBUG /Ox /Zi>
  )
else()

  set(TILEDBSOMA_COMPILE_OPTIONS -Wall -Wextra -D_LIBCPP_TYPEINFO_COMPARISON_IMPLEMENTATION=2)

  if (NOT SPDLOG_LINK_SHARED)
    list(APPEND TILEDBSOMA_COMPILE_OPTIONS -DSPDLOG_USE_STD_FORMAT)
  else()
    add_definitions(-DSPDLOG_COMPILED_LIB)
  endif()

  if(${TILEDBSOMA_ENABLE_WERROR})
    set(TILEDBSOMA_WERROR_OPTION -Werror)
  else()
    set(TILEDBSOMA_WERROR_OPTION "")
  endif()

  # Build-specific flags: Must use generator expressions to support multi-configuration
  # build generators that set the config type at build time.
  list(APPEND
    TILEDBSOMA_COMPILE_OPTIONS
    $<$<CONFIG:Debug>: -DDEBUG -O0 -g3 -ggdb3 -gdwarf-3>
  )
  list(APPEND
    TILEDBSOMA_COMPILE_OPTIONS
    $<$<CONFIG:Release>: -DNDEBUG -O3>
  )
  list(APPEND
    TILEDBSOMA_COMPILE_OPTIONS
    $<$<CONFIG:RelWithDebInfo>: -DNDEBUG -O3 -g3 -ggdb3 -gdwarf-3>
  )
  list(APPEND
    TILEDBSOMA_COMPILE_OPTIONS
    $<$<CONFIG:ASAN,TSAN,LSAN,UBSAN,MSAN>: -DDEBUG -O1 -g -fno-omit-frame-pointer  -fno-optimize-sibling-calls>
  )


  # TODO: There is a bug in CMake 3.22.1 but not 3.27.6 where nested generators
  # are not fully evaluated in the `target_link_library` commands. The GENEX_EVAL
  # commands are not multi-configuration friendly and should be removed after we
  # update CMake.
  # See https://www.reddit.com/r/cmake/comments/17d70h6/why_arent_generator_expressions_evaluated_for/
  # for this placeholder solution.
  set(TILEDBSOMA_SANITIZER_FLAG "")
  list(APPEND TILEDBSOMA_SANITIZER_FLAG "\$<GENEX_EVAL:$<$<CONFIG:ASAN>:-fsanitize=address>>")
  list(APPEND TILEDBSOMA_SANITIZER_FLAG "\$<GENEX_EVAL:$<$<CONFIG:LSAN>:-fsanitize=leak>>")
  list(APPEND TILEDBSOMA_SANITIZER_FLAG "\$<GENEX_EVAL:$<$<CONFIG:TSAN>:-fsanitize=thread>>")
  list(APPEND TILEDBSOMA_SANITIZER_FLAG "\$<GENEX_EVAL:$<$<CONFIG:UBSAN>:-fsanitize=undefined>>")
  list(APPEND TILEDBSOMA_SANITIZER_FLAG "\$<GENEX_EVAL:$<$<CONFIG:MSAN>:-fsanitize=memory>>")


  # Compiler specific additions:
  if(NOT CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    # Use -Wno-literal-suffix on Linux for C++ libtiledbsoma target.
    list(APPEND TILEDBSOMA_COMPILE_OPTIONS -Wno-literal-suffix)
  endif()

endif()


# Definitions for all targets
add_definitions(-D_FILE_OFFSET_BITS=64)

# Disable incorrect availability check on conda build
# https://conda-forge.org/docs/maintainer/knowledge_base.html#newer-c-features-with-old-sdk
add_definitions(-D_LIBCPP_DISABLE_AVAILABILITY)

# AVX2 flag
include(CheckAVX2Support)
CheckAVX2Support()

if(COMPILER_SUPPORTS_AVX2)
  list(APPEND TILEDBSOMA_COMPILE_OPTIONS ${COMPILER_AVX2_FLAG})
endif()


## Leaving this helper. Debugging generator expressions in CMake can be difficult.
## To view compile options uncomment the below and run the target.
## (Hint: This target is inside build/libtiledbsoma)
add_custom_target(
    debugflag COMMAND ${CMAKE_COMMAND} -E echo "compile options: ${TILEDBSOMA_COMPILE_OPTIONS}"
)


# ###########################################################
# Regular build
# ###########################################################

# Adding coverage flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} $ENV{TILEDBSOMA_COVERAGE}")
set(CMAKE_C_FLAGS "${CMAKE_CXX_FLAGS} $ENV{TILEDBSOMA_COVERAGE}")

# Find required dependencies
# See ../cmake/Modules/*.cmake for provenance/version info
find_package(TileDB_EP REQUIRED)
find_package(Spdlog_EP REQUIRED)

add_subdirectory(src)

if(TILEDBSOMA_ENABLE_TESTING)
  enable_testing()
  add_subdirectory(test)
endif()
