#
# CMakeLists.txt for TileDB-SOMA Python API
#
# Licensed under the MIT License.
# Copyright (c) TileDB, Inc. and the Chan Zuckerberg Initiative Foundation
#

cmake_minimum_required(VERSION 3.21)

# Get version from version.py
if(NOT DEFINED SKBUILD_PROJECT_VERSION)
    file(STRINGS "${CMAKE_CURRENT_SOURCE_DIR}/version.py" VERSION_LINE REGEX "^__version__")
    string(REGEX REPLACE "^__version__ = ['\"]([^'\"]+)['\"].*$" "\\1" TILEDBSOMA_VERSION "${VERSION_LINE}")
else()
    set(TILEDBSOMA_VERSION ${SKBUILD_PROJECT_VERSION})
endif()

project(tiledbsoma-python 
    VERSION ${TILEDBSOMA_VERSION}
    LANGUAGES CXX
)

message(STATUS "Building TileDB-SOMA Python API version ${TILEDBSOMA_VERSION}")

# Set C++20 standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Platform-specific settings
if(APPLE)
    set(CMAKE_MACOSX_RPATH ON)
    set(CMAKE_FIND_FRAMEWORK NEVER)
    # Set minimum macOS version
    if(NOT CMAKE_OSX_DEPLOYMENT_TARGET)
        set(CMAKE_OSX_DEPLOYMENT_TARGET "13.3")
    endif()
endif()

# Export all symbols on Windows
set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)

# Find or build libtiledbsoma
# Check if we should build libtiledbsoma or use existing installation
option(BUILD_LIBTILEDBSOMA "Build libtiledbsoma from source" ON)

if(BUILD_LIBTILEDBSOMA)
    message(STATUS "Building libtiledbsoma from source using vcpkg")
    
    set(LIBTILEDBSOMA_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../../libtiledbsoma")
    set(LIBTILEDBSOMA_BUILD_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../../build")
    set(LIBTILEDBSOMA_INSTALL_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../../dist")
    
    # Verify vcpkg is bootstrapped
    set(VCPKG_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/../../vcpkg")
    if(NOT EXISTS "${VCPKG_ROOT}/vcpkg")
        message(FATAL_ERROR "vcpkg not bootstrapped. Please run: ${VCPKG_ROOT}/bootstrap-vcpkg.sh")
    endif()
    
    # Determine vcpkg toolchain file path
    # Use environment variable if set (from setup_build_env.sh), otherwise use local path
    if(DEFINED ENV{CMAKE_TOOLCHAIN_FILE})
        set(VCPKG_TOOLCHAIN_FILE "$ENV{CMAKE_TOOLCHAIN_FILE}")
        message(STATUS "Using vcpkg toolchain from environment: ${VCPKG_TOOLCHAIN_FILE}")
    else()
        set(VCPKG_TOOLCHAIN_FILE "${VCPKG_ROOT}/scripts/buildsystems/vcpkg.cmake")
        message(STATUS "Using vcpkg toolchain from local path: ${VCPKG_TOOLCHAIN_FILE}")
    endif()
    
    if(NOT EXISTS "${VCPKG_TOOLCHAIN_FILE}")
        message(FATAL_ERROR "vcpkg toolchain not found at: ${VCPKG_TOOLCHAIN_FILE}")
    endif()
    
    # Configure libtiledbsoma using vcpkg preset
    # The preset automatically handles vcpkg dependencies (TileDB, spdlog, etc.)
    # Set up environment variables for the vcpkg preset
    set(VCPKG_ENV CMAKE_TOOLCHAIN_FILE=${VCPKG_TOOLCHAIN_FILE})
    
    # Platform-specific environment variables
    if(UNIX AND NOT APPLE AND CMAKE_SYSTEM_PROCESSOR MATCHES "aarch64|arm64")
        # Linux ARM: vcpkg requires using system binaries
        list(APPEND VCPKG_ENV VCPKG_FORCE_SYSTEM_BINARIES=1)
    elseif(DEFINED ENV{VCPKG_FORCE_SYSTEM_BINARIES})
        # Pass through from environment if set (e.g., from cibuildwheel)
        list(APPEND VCPKG_ENV VCPKG_FORCE_SYSTEM_BINARIES=$ENV{VCPKG_FORCE_SYSTEM_BINARIES})
    endif()
    
    if(APPLE AND CMAKE_OSX_DEPLOYMENT_TARGET)
        # macOS: Pass deployment target through
        list(APPEND VCPKG_ENV MACOSX_DEPLOYMENT_TARGET=${CMAKE_OSX_DEPLOYMENT_TARGET})
    elseif(DEFINED ENV{MACOSX_DEPLOYMENT_TARGET})
        # Pass through from environment if set (e.g., from cibuildwheel)
        list(APPEND VCPKG_ENV MACOSX_DEPLOYMENT_TARGET=$ENV{MACOSX_DEPLOYMENT_TARGET})
    endif()
    
    # Build CMake arguments
    set(LIBTILEDBSOMA_CMAKE_ARGS
        --preset vcpkg
        -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
        -DTILEDBSOMA_BUILD_PYTHON=OFF
        -DTILEDBSOMA_BUILD_CLI=OFF
        -DTILEDBSOMA_ENABLE_TESTING=OFF
        -DOVERRIDE_INSTALL_PREFIX=OFF
    )
    
    # Add TileDB deprecation option if specified
    if(DEFINED TILEDB_REMOVE_DEPRECATIONS)
        list(APPEND LIBTILEDBSOMA_CMAKE_ARGS -DTILEDB_REMOVE_DEPRECATIONS=${TILEDB_REMOVE_DEPRECATIONS})
    endif()
    
    # Add OSX-specific arguments if on macOS
    if(APPLE)
        if(CMAKE_OSX_DEPLOYMENT_TARGET)
            list(APPEND LIBTILEDBSOMA_CMAKE_ARGS -DCMAKE_OSX_DEPLOYMENT_TARGET=${CMAKE_OSX_DEPLOYMENT_TARGET})
        endif()
        if(CMAKE_OSX_ARCHITECTURES)
            list(APPEND LIBTILEDBSOMA_CMAKE_ARGS -DCMAKE_OSX_ARCHITECTURES=${CMAKE_OSX_ARCHITECTURES})
        endif()
    endif()
    
    execute_process(
        COMMAND ${CMAKE_COMMAND} -E env ${VCPKG_ENV}
            ${CMAKE_COMMAND} ${LIBTILEDBSOMA_CMAKE_ARGS}
        WORKING_DIRECTORY ${LIBTILEDBSOMA_SOURCE_DIR}
        RESULT_VARIABLE LIBTILEDBSOMA_CONFIGURE_RESULT
        OUTPUT_VARIABLE LIBTILEDBSOMA_CONFIGURE_OUTPUT
        ERROR_VARIABLE LIBTILEDBSOMA_CONFIGURE_OUTPUT
    )
    
    if(NOT LIBTILEDBSOMA_CONFIGURE_RESULT EQUAL 0)
        message(FATAL_ERROR "Failed to configure libtiledbsoma:\n${LIBTILEDBSOMA_CONFIGURE_OUTPUT}")
    endif()
    
    # Build libtiledbsoma
    # vcpkg automatically installs and links all dependencies
    execute_process(
        COMMAND ${CMAKE_COMMAND} --build ${LIBTILEDBSOMA_BUILD_DIR} --config ${CMAKE_BUILD_TYPE} -j
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../..
        RESULT_VARIABLE LIBTILEDBSOMA_BUILD_RESULT
        OUTPUT_VARIABLE LIBTILEDBSOMA_BUILD_OUTPUT
        ERROR_VARIABLE LIBTILEDBSOMA_BUILD_OUTPUT
    )
    
    if(NOT LIBTILEDBSOMA_BUILD_RESULT EQUAL 0)
        message(FATAL_ERROR "Failed to build libtiledbsoma:\n${LIBTILEDBSOMA_BUILD_OUTPUT}")
    endif()
    
    # Install libtiledbsoma headers and libraries
    execute_process(
        COMMAND ${CMAKE_COMMAND} --build ${LIBTILEDBSOMA_BUILD_DIR} --target install-libtiledbsoma --config ${CMAKE_BUILD_TYPE}
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../..
        RESULT_VARIABLE LIBTILEDBSOMA_INSTALL_RESULT
        OUTPUT_VARIABLE LIBTILEDBSOMA_INSTALL_OUTPUT
        ERROR_VARIABLE LIBTILEDBSOMA_INSTALL_OUTPUT
    )
    
    if(NOT LIBTILEDBSOMA_INSTALL_RESULT EQUAL 0)
        message(FATAL_ERROR "Failed to install libtiledbsoma:\n${LIBTILEDBSOMA_INSTALL_OUTPUT}")
    endif()
    
    # Point to the installed libtiledbsoma
    set(LIBTILEDBSOMA_INCLUDE_DIR "${LIBTILEDBSOMA_INSTALL_DIR}/include")
    set(LIBTILEDBSOMA_LIB_DIR "${LIBTILEDBSOMA_INSTALL_DIR}/lib")
    
    # External dependencies (thread_pool, etc.) are in the libtiledbsoma source tree
    set(LIBTILEDBSOMA_EXTERNALS_INCLUDE "${LIBTILEDBSOMA_SOURCE_DIR}/src/external/include")
    
    # Find vcpkg triplet directory for libraries
    file(GLOB VCPKG_TRIPLET_DIRS "${LIBTILEDBSOMA_BUILD_DIR}/vcpkg_installed/*")
    foreach(TRIPLET_DIR ${VCPKG_TRIPLET_DIRS})
        if(IS_DIRECTORY "${TRIPLET_DIR}" AND EXISTS "${TRIPLET_DIR}/lib")
            get_filename_component(TRIPLET_NAME "${TRIPLET_DIR}" NAME)
            if(NOT TRIPLET_NAME STREQUAL "vcpkg")
                set(LIBTILEDBSOMA_EXTERNALS_LIB "${TRIPLET_DIR}/lib")
                break()
            endif()
        endif()
    endforeach()
    
else()
    # Use existing libtiledbsoma installation
    find_package(TileDBSOMA REQUIRED)
    set(LIBTILEDBSOMA_INSTALL_DIR "${TileDBSOMA_DIR}/..")
    set(LIBTILEDBSOMA_INCLUDE_DIR "${TileDBSOMA_INCLUDE_DIR}")
    set(LIBTILEDBSOMA_LIB_DIR "${LIBTILEDBSOMA_INSTALL_DIR}/lib")
endif()

# Find Python and pybind11
find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)
find_package(pybind11 CONFIG REQUIRED)

# Compiler flags
if(MSVC)
    set(CXX_FLAGS /W4 /wd4101 /wd4146 /wd4244 /wd4251 /wd4456 /wd4457 /wd4702 /wd4800 /wd4996 /DNOGDI)
else()
    set(CXX_FLAGS -Wall -Wextra -D_LIBCPP_TYPEINFO_COMPARISON_IMPLEMENTATION=2)
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
        list(APPEND CXX_FLAGS -Wno-literal-suffix)
    endif()
endif()

# Platform-specific compiler flags
if(APPLE)
    list(APPEND CXX_FLAGS -mmacosx-version-min=13.3)
endif()

# AVX2 support
include(CheckCXXCompilerFlag)
check_cxx_compiler_flag("-mavx2" COMPILER_SUPPORTS_AVX2)
if(COMPILER_SUPPORTS_AVX2 AND CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64|AMD64")
    list(APPEND CXX_FLAGS -mavx2)
endif()

# Create the Python extension module
pybind11_add_module(pytiledbsoma
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

# Configure include directories for Python extension
# Include libtiledbsoma headers and external dependencies (thread_pool, etc.)
target_include_directories(pytiledbsoma PRIVATE
    ${LIBTILEDBSOMA_INCLUDE_DIR}              # Main libtiledbsoma headers
    ${LIBTILEDBSOMA_INCLUDE_DIR}/tiledbsoma   # Subdirectory headers
)

# Add external dependencies include path (thread_pool, nanoarrow, etc.)
# These are vendored in the libtiledbsoma source tree
if(BUILD_LIBTILEDBSOMA)
    target_include_directories(pytiledbsoma PRIVATE ${LIBTILEDBSOMA_EXTERNALS_INCLUDE})
endif()

# Link Python extension against libtiledbsoma
target_link_directories(pytiledbsoma PRIVATE ${LIBTILEDBSOMA_LIB_DIR})
target_link_libraries(pytiledbsoma PRIVATE tiledbsoma)

# Set compiler options and C++20 standard
target_compile_options(pytiledbsoma PRIVATE ${CXX_FLAGS})
target_compile_features(pytiledbsoma PRIVATE cxx_std_20)

# Configure RPATH for runtime library loading
# This ensures the Python extension can find shared libraries bundled in the wheel
if(APPLE)
    # macOS: Look in same directory and .dylibs subdirectory
    set_target_properties(pytiledbsoma PROPERTIES
        BUILD_WITH_INSTALL_RPATH TRUE
        INSTALL_RPATH "@loader_path;@loader_path/../.dylibs"
    )
elseif(UNIX)
    # Linux: Look in same directory and .libs subdirectory
    set_target_properties(pytiledbsoma PROPERTIES
        BUILD_WITH_INSTALL_RPATH TRUE
        INSTALL_RPATH "$ORIGIN:$ORIGIN/../.libs"
    )
endif()

# Install the extension module
install(TARGETS pytiledbsoma
    LIBRARY DESTINATION tiledbsoma
    RUNTIME DESTINATION tiledbsoma
)

# Install libtiledbsoma and vcpkg dependencies alongside the Python module
# This bundles all required shared libraries into the wheel for distribution
if(BUILD_LIBTILEDBSOMA)
    if(APPLE)
        # Install libtiledbsoma dylib
        install(DIRECTORY ${LIBTILEDBSOMA_LIB_DIR}/
            DESTINATION tiledbsoma/.dylibs
            FILES_MATCHING 
            PATTERN "*.dylib"
            PATTERN "*.dylib.*"
            PATTERN "cmake" EXCLUDE
            PATTERN "pkgconfig" EXCLUDE
        )
        # Install vcpkg-provided dependencies (TileDB, spdlog, etc.)
        if(EXISTS "${LIBTILEDBSOMA_EXTERNALS_LIB}")
            install(DIRECTORY ${LIBTILEDBSOMA_EXTERNALS_LIB}/
                DESTINATION tiledbsoma/.dylibs
                FILES_MATCHING 
                PATTERN "*.dylib"
                PATTERN "*.dylib.*"
                PATTERN "cmake" EXCLUDE
                PATTERN "pkgconfig" EXCLUDE
            )
        endif()
    elseif(UNIX)
        # Install libtiledbsoma shared object
        install(DIRECTORY ${LIBTILEDBSOMA_LIB_DIR}/
            DESTINATION tiledbsoma/.libs
            FILES_MATCHING 
            PATTERN "*.so"
            PATTERN "*.so.*"
            PATTERN "cmake" EXCLUDE
            PATTERN "pkgconfig" EXCLUDE
        )
        # Install vcpkg-provided dependencies (TileDB, spdlog, etc.)
        if(EXISTS "${LIBTILEDBSOMA_EXTERNALS_LIB}")
            install(DIRECTORY ${LIBTILEDBSOMA_EXTERNALS_LIB}/
                DESTINATION tiledbsoma/.libs
                FILES_MATCHING 
                PATTERN "*.so"
                PATTERN "*.so.*"
                PATTERN "cmake" EXCLUDE
                PATTERN "pkgconfig" EXCLUDE
            )
        endif()
    elseif(WIN32)
        # Install libtiledbsoma DLL
        install(DIRECTORY ${LIBTILEDBSOMA_INSTALL_DIR}/bin/
            DESTINATION tiledbsoma
            FILES_MATCHING 
            PATTERN "*.dll"
        )
        # Install vcpkg-provided dependencies (TileDB, spdlog, etc.)
        if(EXISTS "${LIBTILEDBSOMA_EXTERNALS_LIB}/../bin")
            install(DIRECTORY ${LIBTILEDBSOMA_EXTERNALS_LIB}/../bin/
                DESTINATION tiledbsoma
                FILES_MATCHING 
                PATTERN "*.dll"
            )
        endif()
    endif()
endif()

# Print configuration summary
message(STATUS "===================================")
message(STATUS "TileDB-SOMA Python Build Configuration:")
message(STATUS "  Version: ${TILEDBSOMA_VERSION}")
message(STATUS "  Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "  C++ compiler: ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")
message(STATUS "  Python: ${Python_VERSION}")
message(STATUS "  Build libtiledbsoma: ${BUILD_LIBTILEDBSOMA}")
message(STATUS "  libtiledbsoma include: ${LIBTILEDBSOMA_INCLUDE_DIR}")
message(STATUS "  libtiledbsoma install: ${LIBTILEDBSOMA_INSTALL_DIR}")
if(BUILD_LIBTILEDBSOMA)
    message(STATUS "  vcpkg externals include: ${LIBTILEDBSOMA_EXTERNALS_INCLUDE}")
    message(STATUS "  vcpkg externals lib: ${LIBTILEDBSOMA_EXTERNALS_LIB}")
endif()
if(DEFINED CMAKE_TOOLCHAIN_FILE)
    message(STATUS "  vcpkg toolchain: ${CMAKE_TOOLCHAIN_FILE}")
endif()
message(STATUS "===================================")

