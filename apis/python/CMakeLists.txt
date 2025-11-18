cmake_minimum_required(VERSION 3.21)
cmake_policy(SET CMP0148 NEW) 

# Set macOS architecture before project() to avoid vcpkg warnings
if(APPLE AND NOT DEFINED CMAKE_OSX_ARCHITECTURES)
    # Detect system architecture
    execute_process(
        COMMAND uname -m
        OUTPUT_VARIABLE DETECTED_ARCH
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    set(CMAKE_OSX_ARCHITECTURES "${DETECTED_ARCH}" CACHE STRING "macOS architecture")
    message(STATUS "Auto-detected macOS architecture: ${CMAKE_OSX_ARCHITECTURES}")
endif()

project(tiledbsoma-python LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find pybind11 - required for pybind11_add_module
find_package(pybind11 REQUIRED)
find_package(Ninja)

# Find Python
find_package(Python COMPONENTS Interpreter Development REQUIRED)

# Get version from setuptools_scm generated file or fallback
# setuptools_scm writes version to src/tiledbsoma/_version.py if configured
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/tiledbsoma/_version.py")
    # Try to read version from setuptools_scm generated file
    file(READ "${CMAKE_CURRENT_SOURCE_DIR}/tiledbsoma/_version.py" VERSION_PY_CONTENT)
    string(REGEX MATCH "__version__ = \"([^\"]+)\"" _ "${VERSION_PY_CONTENT}")
    if(CMAKE_MATCH_1)
        set(VERSION_STRING "${CMAKE_MATCH_1}")
        message(STATUS "Package version (from setuptools_scm): ${VERSION_STRING}")
    else()
        set(VERSION_STRING "0.0.0.dev0")
        message(WARNING "Could not parse version from _version.py, using default")
    endif()
elseif(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/RELEASE-VERSION")
    # Fallback to RELEASE-VERSION file if it exists (matches version.py convention)
    file(READ "${CMAKE_CURRENT_SOURCE_DIR}/RELEASE-VERSION" VERSION_STRING)
    string(STRIP "${VERSION_STRING}" VERSION_STRING)
    message(STATUS "Package version (from RELEASE-VERSION file): ${VERSION_STRING}")
else()
    set(VERSION_STRING "0.0.0.dev0")
    message(WARNING "No version file found, using default version")
endif()

# Build libtiledbsoma using ExternalProject
# This will build libtiledbsoma and all its dependencies via vcpkg
# Check for the sdist-vendored directory first
set(LIBTDB_SDIST_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/dist_links/libtiledbsoma")

if(EXISTS "${LIBTDB_SDIST_SOURCE_DIR}")
    # We are in an sdist build. Use the vendored path.
    message(STATUS "Using vendored libtiledbsoma source from sdist")
    set(LIBTILEDBSOMA_SOURCE_DIR "${LIBTDB_SDIST_SOURCE_DIR}")
else()
    # We are in a dev (Git) build. Use the relative path.
    message(STATUS "Using libtiledbsoma source from Git repo")
    set(LIBTILEDBSOMA_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../../libtiledbsoma")
endif()
set(LIBTILEDBSOMA_BINARY_DIR "${CMAKE_CURRENT_BINARY_DIR}/libtiledbsoma-build")
set(LIBTILEDBSOMA_INSTALL_PREFIX "${CMAKE_CURRENT_BINARY_DIR}/install")

# Include ExternalProject
include(ExternalProject)

# Prepare CMake configure arguments for libtiledbsoma
set(LIBTILEDBSOMA_CMAKE_ARGS
    -DCMAKE_BUILD_TYPE=Release
    -DCMAKE_INSTALL_PREFIX=${LIBTILEDBSOMA_INSTALL_PREFIX}
    -DOVERRIDE_INSTALL_PREFIX=OFF
    -DTILEDB_REMOVE_DEPRECATIONS=ON
    -DSUPERBUILD=ON
)

# Forward compilers if explicitly set
if(CMAKE_C_COMPILER)
    list(APPEND LIBTILEDBSOMA_CMAKE_ARGS
        -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
    )
endif()

if(CMAKE_CXX_COMPILER)
    list(APPEND LIBTILEDBSOMA_CMAKE_ARGS
        -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
    )
endif()

if(NINJA_FOUND)
    list(APPEND LIBTILEDBSOMA_CMAKE_ARGS
        -DCMAKE_MAKE_PROGRAM=${NINJA_EXECUTABLE}
    )
endif()

# Forward CMAKE_TOOLCHAIN_FILE - critical for vcpkg
message(STATUS "CMAKE_TOOLCHAIN_FILE: ${CMAKE_TOOLCHAIN_FILE}")
message(STATUS "ENV CMAKE_TOOLCHAIN_FILE: $ENV{CMAKE_TOOLCHAIN_FILE}")
if(DEFINED CMAKE_TOOLCHAIN_FILE AND CMAKE_TOOLCHAIN_FILE)
    message(STATUS "Adding CMAKE_TOOLCHAIN_FILE from CMake variable: ${CMAKE_TOOLCHAIN_FILE}")
    list(APPEND LIBTILEDBSOMA_CMAKE_ARGS
        -DCMAKE_TOOLCHAIN_FILE=${CMAKE_TOOLCHAIN_FILE}
    )
elseif(DEFINED ENV{CMAKE_TOOLCHAIN_FILE} AND NOT "$ENV{CMAKE_TOOLCHAIN_FILE}" STREQUAL "")
    message(STATUS "Adding CMAKE_TOOLCHAIN_FILE from environment: $ENV{CMAKE_TOOLCHAIN_FILE}")
    list(APPEND LIBTILEDBSOMA_CMAKE_ARGS
        -DCMAKE_TOOLCHAIN_FILE=$ENV{CMAKE_TOOLCHAIN_FILE}
    )
else()
    message(WARNING "CMAKE_TOOLCHAIN_FILE not set! vcpkg will not be used.")
endif()

# Forward macOS architecture if set
if(DEFINED CMAKE_OSX_ARCHITECTURES)
    list(APPEND LIBTILEDBSOMA_CMAKE_ARGS
        -DCMAKE_OSX_ARCHITECTURES=${CMAKE_OSX_ARCHITECTURES}
    )
endif()

# Forward vcpkg triplet if set (from environment or CMake variable)
if(VCPKG_TARGET_TRIPLET)
    list(APPEND LIBTILEDBSOMA_CMAKE_ARGS
        -DVCPKG_TARGET_TRIPLET=${VCPKG_TARGET_TRIPLET}
    )
elseif(DEFINED ENV{VCPKG_TARGET_TRIPLET})
    list(APPEND LIBTILEDBSOMA_CMAKE_ARGS
        -DVCPKG_TARGET_TRIPLET=$ENV{VCPKG_TARGET_TRIPLET}
    )
endif()

# Forward VCPKG_FORCE_SYSTEM_BINARIES for ARM builds
if(DEFINED ENV{VCPKG_FORCE_SYSTEM_BINARIES})
    list(APPEND LIBTILEDBSOMA_CMAKE_ARGS
        -DVCPKG_FORCE_SYSTEM_BINARIES=$ENV{VCPKG_FORCE_SYSTEM_BINARIES}
    )
endif()

# Prepare environment for ExternalProject (to pass vcpkg environment)
message(STATUS "ENV VCPKG_ROOT: $ENV{VCPKG_ROOT}")
message(STATUS "ENV VCPKG_FORCE_SYSTEM_BINARIES: $ENV{VCPKG_FORCE_SYSTEM_BINARIES}")
set(LIBTILEDBSOMA_ENV_VARS "")
if(DEFINED ENV{VCPKG_ROOT} AND NOT "$ENV{VCPKG_ROOT}" STREQUAL "")
    list(APPEND LIBTILEDBSOMA_ENV_VARS "VCPKG_ROOT=$ENV{VCPKG_ROOT}")
    message(STATUS "Adding VCPKG_ROOT to environment: $ENV{VCPKG_ROOT}")
endif()
if(DEFINED ENV{VCPKG_FORCE_SYSTEM_BINARIES} AND NOT "$ENV{VCPKG_FORCE_SYSTEM_BINARIES}" STREQUAL "")
    list(APPEND LIBTILEDBSOMA_ENV_VARS "VCPKG_FORCE_SYSTEM_BINARIES=$ENV{VCPKG_FORCE_SYSTEM_BINARIES}")
    message(STATUS "Adding VCPKG_FORCE_SYSTEM_BINARIES to environment: $ENV{VCPKG_FORCE_SYSTEM_BINARIES}")
endif()
message(STATUS "Environment variables for ExternalProject: ${LIBTILEDBSOMA_ENV_VARS}")

# Build libtiledbsoma as ExternalProject
# Note: We configure libtiledbsoma with SUPERBUILD=ON so it builds TileDB via superbuild
# while using vcpkg for other dependencies (spdlog, fmt, catch2)
# Use CONFIGURE_COMMAND to pass environment variables during configuration
ExternalProject_Add(
    libtiledbsoma_ep
    SOURCE_DIR ${LIBTILEDBSOMA_SOURCE_DIR}
    BINARY_DIR ${LIBTILEDBSOMA_BINARY_DIR}
    CONFIGURE_COMMAND ${CMAKE_COMMAND} -E env ${LIBTILEDBSOMA_ENV_VARS}
        ${CMAKE_COMMAND} -G ${CMAKE_GENERATOR}
        ${LIBTILEDBSOMA_CMAKE_ARGS}
        -S ${LIBTILEDBSOMA_SOURCE_DIR}
        -B ${LIBTILEDBSOMA_BINARY_DIR}
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
    BUILD_ALWAYS OFF
)

# Custom command to build and install libtiledbsoma after configuration
# Pass environment variables for vcpkg
ExternalProject_Add_Step(libtiledbsoma_ep build_and_install
    COMMAND ${CMAKE_COMMAND} -E env ${LIBTILEDBSOMA_ENV_VARS} ${CMAKE_COMMAND} --build ${LIBTILEDBSOMA_BINARY_DIR}
    COMMAND ${CMAKE_COMMAND} -E env ${LIBTILEDBSOMA_ENV_VARS} ${CMAKE_COMMAND} --build ${LIBTILEDBSOMA_BINARY_DIR}/libtiledbsoma --target install
    DEPENDEES configure
    DEPENDERS build
    ALWAYS OFF
)

# Note: We can't use find_package here because libtiledbsoma_ep hasn't built yet
# Instead, we'll set up the paths directly and link manually
# The ExternalProject will build libtiledbsoma before pytiledbsoma (via add_dependencies)

# Include directories
set(INCLUDE_DIRS
    ${LIBTILEDBSOMA_INSTALL_PREFIX}/include
    ${LIBTILEDBSOMA_INSTALL_PREFIX}/include/tiledbsoma
    ${CMAKE_CURRENT_SOURCE_DIR}/src/tiledbsoma
)

# Library directories
set(LIB_DIRS
    ${LIBTILEDBSOMA_INSTALL_PREFIX}/lib
    ${LIBTILEDBSOMA_INSTALL_PREFIX}/lib64
)

# Source files for Python extension
set(PYTHON_SOURCES
    src/tiledbsoma/common.cc
    src/tiledbsoma/fastercsx.cc
    src/tiledbsoma/reindexer.cc
    src/tiledbsoma/query_condition.cc
    src/tiledbsoma/soma_vfs.cc
    src/tiledbsoma/soma_context.cc
    src/tiledbsoma/soma_array.cc
    src/tiledbsoma/soma_object.cc
    src/tiledbsoma/soma_column.cc
    src/tiledbsoma/soma_dataframe.cc
    src/tiledbsoma/soma_point_cloud_dataframe.cc
    src/tiledbsoma/soma_geometry_dataframe.cc
    src/tiledbsoma/soma_dense_ndarray.cc
    src/tiledbsoma/soma_sparse_ndarray.cc
    src/tiledbsoma/soma_group.cc
    src/tiledbsoma/soma_collection.cc
    src/tiledbsoma/coordinate_selection.cc
    src/tiledbsoma/managed_query.cc
    src/tiledbsoma/transformer.cc
    src/tiledbsoma/pytiledbsoma.cc
)

# Compiler flags
set(CXX_FLAGS "-O3")

# AVX2 support for x86_64
if(CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64" OR CMAKE_SYSTEM_PROCESSOR MATCHES "AMD64")
    list(APPEND CXX_FLAGS "-mavx2")
endif()

# macOS-specific flags
if(APPLE)
    list(APPEND CXX_FLAGS "-mmacosx-version-min=13.0")
    list(APPEND CXX_FLAGS "-D_LIBCPP_TYPEINFO_COMPARISON_IMPLEMENTATION=2")
endif()

# Create Python extension module
pybind11_add_module(
    pytiledbsoma
    ${PYTHON_SOURCES}
)

# Make sure libtiledbsoma is built before we try to link to it
add_dependencies(pytiledbsoma libtiledbsoma_ep)

target_include_directories(pytiledbsoma PRIVATE ${INCLUDE_DIRS})
target_link_directories(pytiledbsoma PRIVATE ${LIB_DIRS})
target_link_libraries(pytiledbsoma PRIVATE tiledbsoma)

# Set compiler flags
target_compile_options(pytiledbsoma PRIVATE ${CXX_FLAGS})

# Set C++ standard and RPATH for finding bundled libraries
if(APPLE)
    # macOS uses @loader_path instead of $ORIGIN
    set(RPATH_ORIGIN "@loader_path")
else()
    # Linux uses $ORIGIN
    set(RPATH_ORIGIN "$ORIGIN")
endif()

set_target_properties(pytiledbsoma PROPERTIES
    CXX_STANDARD 20
    CXX_STANDARD_REQUIRED ON
    # Set RPATH to look in same directory as extension module
    BUILD_RPATH "${LIBTILEDBSOMA_INSTALL_PREFIX}/lib;${LIBTILEDBSOMA_INSTALL_PREFIX}/lib64"
    INSTALL_RPATH "${RPATH_ORIGIN}"
    INSTALL_RPATH_USE_LINK_PATH FALSE
    BUILD_WITH_INSTALL_RPATH FALSE
    SKIP_BUILD_RPATH FALSE
)

# Install the extension module
install(TARGETS pytiledbsoma DESTINATION tiledbsoma)

# Debug: Check if libraries exist before installing
if(APPLE)
    message(STATUS "Checking for: ${LIBTILEDBSOMA_INSTALL_PREFIX}/lib/libtiledbsoma.dylib")
elseif(UNIX)
    message(STATUS "Checking for: ${LIBTILEDBSOMA_INSTALL_PREFIX}/lib/libtiledbsoma.so")
    message(STATUS "Checking for: ${LIBTILEDBSOMA_INSTALL_PREFIX}/lib64/libtiledbsoma.so")
endif()

# Install shared libraries - try both lib and lib64 locations
if(APPLE)
    install(FILES ${LIBTILEDBSOMA_INSTALL_PREFIX}/lib/libtiledbsoma.dylib
        DESTINATION tiledbsoma
    )
elseif(UNIX)
    # Install from whichever location exists (lib or lib64)
    install(CODE "
        if(EXISTS \"${LIBTILEDBSOMA_INSTALL_PREFIX}/lib/libtiledbsoma.so\")
            file(INSTALL \"${LIBTILEDBSOMA_INSTALL_PREFIX}/lib/libtiledbsoma.so\" 
                 DESTINATION \"\${CMAKE_INSTALL_PREFIX}/tiledbsoma\")
            message(STATUS \"Installed libtiledbsoma.so from lib/\")
        elseif(EXISTS \"${LIBTILEDBSOMA_INSTALL_PREFIX}/lib64/libtiledbsoma.so\")
            file(INSTALL \"${LIBTILEDBSOMA_INSTALL_PREFIX}/lib64/libtiledbsoma.so\" 
                 DESTINATION \"\${CMAKE_INSTALL_PREFIX}/tiledbsoma\")
            message(STATUS \"Installed libtiledbsoma.so from lib64/\")
        else()
            message(FATAL_ERROR \"Could not find libtiledbsoma.so in ${LIBTILEDBSOMA_INSTALL_PREFIX}/lib or ${LIBTILEDBSOMA_INSTALL_PREFIX}/lib64\")
        endif()
    ")
elseif(WIN32)
    install(FILES ${LIBTILEDBSOMA_INSTALL_PREFIX}/bin/tiledbsoma.dll
        DESTINATION tiledbsoma
        OPTIONAL
    )
endif()

# Install TileDB shared libraries if not provided externally
if(NOT DEFINED ENV{TILEDB_PATH})
    if(APPLE)
        install(FILES ${LIBTILEDBSOMA_INSTALL_PREFIX}/lib/libtiledb.dylib
            DESTINATION tiledbsoma
        )
    elseif(UNIX)
        install(CODE "
            if(EXISTS \"${LIBTILEDBSOMA_INSTALL_PREFIX}/lib/libtiledb.so\")
                file(INSTALL \"${LIBTILEDBSOMA_INSTALL_PREFIX}/lib/libtiledb.so\" 
                     DESTINATION \"\${CMAKE_INSTALL_PREFIX}/tiledbsoma\")
                message(STATUS \"Installed libtiledb.so from lib/\")
            elseif(EXISTS \"${LIBTILEDBSOMA_INSTALL_PREFIX}/lib64/libtiledb.so\")
                file(INSTALL \"${LIBTILEDBSOMA_INSTALL_PREFIX}/lib64/libtiledb.so\" 
                     DESTINATION \"\${CMAKE_INSTALL_PREFIX}/tiledbsoma\")
                message(STATUS \"Installed libtiledb.so from lib64/\")
            else()
                message(FATAL_ERROR \"Could not find libtiledb.so in ${LIBTILEDBSOMA_INSTALL_PREFIX}/lib or ${LIBTILEDBSOMA_INSTALL_PREFIX}/lib64\")
            endif()
        ")
    elseif(WIN32)
        install(FILES ${LIBTILEDBSOMA_INSTALL_PREFIX}/bin/tiledb.dll
            DESTINATION tiledbsoma
        )
    endif()
endif()

